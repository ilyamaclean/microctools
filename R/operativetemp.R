#' Calculates angle between organism and sun
#' @description Calculates angle between organism's longitudinal axis and the direction
#' of the solar beam
#' @param alt solar altitude (degrees) as returned by [solalt()]
#' @param azi solar azimuth (degrees) as returned by [solazi]
#' @param slope slope angle (degrees) of organism's longitudinal axis (see details)
#' @param aspect aspect (degrees relative to north) of an organism's longitudinal axis (see details)
#' @param shape logical indicating whether organism is flat (leaf) or not (see details)
#' @return angle between organism and sun (degrees)
#' @seealso [Fprad()], [Radabs()]
#' @export
#' @details
#' Slope and aspect values are relative to the horizontal and north respectively. E.g. an upright
#' human has a slope of 90 degrees, whereas a caterpillar facing down-slope as the same slope
#' and apsect as the underlying terrain. Shape must be specified, as for a flat plate, `azi`
#' has no effect on the angle between the organism and sun when `slope` is 0, but for spheres
#' and cylinders `azi` has no effect on the angle between organism and sun when `slope` is 90.
#' 90 degrees for a human, 0 for an organism on a flat surface. The function `thetangle` is used
#' in the calculation of view factors for direct beam radiation and absorbed radiation and if the orientation of the
#' organism is unknown, absorbed radiation can be calculated by assuming  no direction bias
#' in orientation and deriving results numerically as in the example for [Fprad()].
#' @examples
#' alt <- c(0:90) # solar altitude
#' # Angle for upright human
#' theta <- thetangle(alt, azi = 0, slope = 90, aspect = 0)
#' plot(theta ~ alt, type = "l")
thetangle <- function(alt, azi, slope = 0, aspect = 180, plate = FALSE) {
  if (plate) {
    zen <- (pi/2)-alt*(pi/180)
    azi <- azi*(pi/180)
    sl<-slope*(pi/180)
    aspe <- aspect*(pi/180)
    index <- cos(zen)*cos(sl)+sin(zen)*sin(sl)*cos(azi-aspe)
    theta<-acos(index)*(180/pi)
  } else {
    alt <- (90-alt)*(pi/180)
    azi <- azi*(pi/180)
    sl<-(pi/2)-slope*(pi/180)
    aspe <- aspect*(pi/180)
    index <- cos(alt)*cos(sl)+sin(alt)*sin(sl)*cos(azi-aspe)
    theta<-acos(index)*(180/pi)
  }
  theta
}
#' Calculates view factor of organism for beam radiation
#' @description the view factor is the average flux density of beam radiation over the entire
#' surface of the object relative to the flux density on a flat absobing surface facing the
#' source
#' @param theta angle between the longitudinal axis of the animal and solar beam (degrees) as
#' returned by [thetangle()]
#' @param l length of organism (m) along longitudinal axis
#' @param d diameter of organism (m)
#' @param shape one of `spheroid` (prolate spheroid), `cylinder1` (cylinder with flat ends) or
#' `cylinder2` (cylinder with hemisphere ends) or `plate` (flat plate - e.g. leaf)
#' @return The view factor (ranges between zero and 1, typically ~0.25)
#' @seealso [thetangle()], [Radabs()]
#' @details The view factor is used in calculation of absorbed radiation. It is effectively
#' the ratio of of the projected area in the direction of the solar beam to the total surface
#' area of the object. For a flat plate normal to the beam, the view factor is 0.5, as only have the
#' plate faces the beam.
#' @export
#' @examples
#' # caterpillar on slope (orientation unknown, so calculated by averaging across all orientations)
#' aspect <- c(0:359)
#' theta<-thetangle(alt = 50, azi = 180, slope = 20, aspect)
#' Fp <- Fprad(theta, 0.05, 0.01, shape = "cylinder2")
#' plot(Fp~aspect,type="l")
#' Fpm <- mean(Fp) # mean viewing angle
#' Fpm
Fprad <- function(theta, l, d, shape = "spheroid") {
  if (shape == "spheroid") {
    if (d > l) stop ("for spheroids l must be greater than d")
    x <- d / l
    if (x < 1) {
      btm <- 2*x+((2 * asin(sqrt(1-x^2)))/(sqrt(1-x^2)))
      top <- sqrt(1+(x^2-1)*cos(theta*(pi/180)))
      Fp <- top/btm
    } else Fp <- rep(0.25, length(theta))
  }
  if (shape == "cylinder1") {
    btm <- 2 + (4 * l) / d
    top <- cos(theta*(pi/180))+(4*l*sin(theta*(pi/180)))/(pi*d)
    Fp <- top/btm
  }
  if (shape == "cylinder2") {
    btm <- 4 + (4 * l) / d
    top <- 1 + (4*l*sin(theta*(pi/180)))/(pi*d)
    Fp <- top/btm
  }
  if (shape == "plate") {
    Fp <- cos(theta * (pi/180))
  }
  Fp
}
#' Calculates radiation absorbed by an organism
#' @description calculates the sum of long and shortwave radiation absorbed by and organism
#' @param dni Flux density of beam radiation perpendicular to the solar beam (W / m2)
#' @param dif Flux density of beam radiation perpendicular to the solar beam (W / m2) (see [difprop()])
#' @param si Solar coefficient as returned by [solarcoef()]
#' @param tair temperature of air surrounding organism (deg C)
#' @param tground temperature of ground surface  (deg C)
#' @param l length or organism (m) along longitudinal axis
#' @param d diameter (m) of organism
#' @param theta angle between longitudinal axis of organism and sun as returned by [thetangle]
#' @param refsw reflectivity of organism (shortwave radiation)
#' @param em emissivity of organism (longwave radiation)
#' @param lwta fraction of longwave radiation transmitted downward (see [canlw()].
#' 1 - sky emissivity if no canopy cover.
#' @param lwtg fraction of longwave radiation transmitted upward (see [canlw()]. Determined by canopy
#' below organism. Ignored if `suspended` = FALSE.
#' @param shape  one of `spheroid` (prolate spheroid), `cylinder1` (cylinder with flat ends) or
#' `cylinder2` (cylinder with hemisphere ends) or `plate` (flat plate - e.g. leaf)
#' @param suspended optional logical indicating whether organism is suspended in free air (single
#' leaves, organisms with spindly legs) or on the surface (whole canopies, organisms resting on the ground)
#' @return Flux density of absorbed radiation (see details)
#' @details The returned value is the flux density (i.e. radiation absorbed per metre squared).
#' To convert to total flux, multiply by surface area of organism.
#' @export
Radabs <- function(dni, dif, si, tair, tground, l, d, theta, refsw = 0.5, em = 0.97, lwta = 0.5, lwtg = 1, galb = 0.15,
                   shape = "spheroid", suspended = FALSE) {
  Fp <- Fprad(theta, l, d, shape)
  Lwa <- lwta * 5.67*10^-8*0.97*(tair + 273.15)^4
  Lwg <- lwtg * 5.67*10^-8*0.97*(tground + 273.15)^4
  # View factors: Fp = beam; Fd = diffuse; Fr = reflected; Fa = atmosphere (lw); Fg = ground (lw)
  Sr <- (1-galb)*(si*dni+dif)
  if (suspended) {
    Rabs <- (1-refsw)*(Fp*dni+0.5*dif+0.5*Sr)+em*(0.5*Lwa+0.5*Lwg)
  } else  Rabs <- (1-refsw)*(Fp*dni+dif)+em*Lwa
  Rabs
}
#' Calculates Operative Temperature
#' @description The function `OperativeT` calculates the humid operative temperature of
#' an organism
#' @param Rabs Absorbed long- and shortwave radiation as returned by [Radbs()] (W / m^2)
#' @param gH conductance for heat (mol / m^2 / sec)
#' @param gV conductance for vapour (mol / m^2 / sec)
#' @param tair temperature of air surrounding organism (deg C)
#' @param ea Vapour pressur eof air surrounding organism (kPa)
#' @param pk Atmospheric pressure (kPa)
#' @param surfwet Proportion of surface area acting like a free‐water surface
#' @param em thermal emissivity
#' @param fluxother Fluxe sother than radiation, sensible and latent heat operating on
#' organism (W / m^2). Positive if warming organism (e.g. metabolism), negative if leaving
#' organism (e.g. ground heat flux)
#' @param method one of `iter`, `Penman` or `PenMon` (see details)
#' @return Humid Operative Temperature of organism - i.e. steady-state temperature of
#' organism (deg C)
#' @export
#' @details The Humid Operative Temperature is the temperature at which all fluxes operating
#' on the organism sum to zero. If method = `iter`, this temperature is calaculated iteratively
#' and the inputs must be single values rather than vectors. If method = `open`, radiation
#' emmited is linearized by computing radiative conductance, and the Penman linearization for
#' latent heat flux is used. The results obtained are similar, though not quite as accurate
#' as those obtained using `iter`.If method = `PenmMon`, the Penman-Monteith equation is
#' used, i.e. it is assumed that the temperature and humidity of the environment are influenced
#' by evapotransporation from the organism, such that humidity increases and the air cools.
#' In consequence, there is less latent heat flux, and the operative temperatures are higher.
#' @examples
#' Rabs <- c(200:700)
#' OT1 <- 0
#' for (i in 1:length(Rabs)) OT1[i] <- OperativeT(Rabs[i], 0.3, 0.2, 15, 1.2, method = "iter")
#' OT2 <- OperativeT(Rabs, 0.3, 0.2, 15, 1.2, method = "open")
#' OT3 <- OperativeT(Rabs, 0.3, 0.2, 15, 1.2, method = "PenmMon")
#' plot(OT1 ~ Rabs, type = "l", ylim = c(0,40), lwd = 2)
#' par(new=T)
#' plot(OT2 ~ Rabs, type = "l", ylim = c(0,40), col = "blue", lwd = 2)
#' par(new=T)
#' plot(OT3 ~ Rabs, type = "l", ylim = c(0,40), col = "red", lwd = 2)
OperativeT <- function(Rabs, gH, gV, tair, ea, pk = 101.3, surfwet = 1, em = 0.97,
                       fluxother = 0, method = "iter") {
  Balance <- function(Tb, sb, em, cp, gH, gV, lambda, surfwet, ea, pk, fluxother) {
    estb <- 0.6108*exp(17.27*Tb/(Tb+237.3))
    Rem <- sb*em*(Tb+273.15)^4
    HH <- cp*gH*(Tb-tair)
    LL <- surfwet*(lambda/pk)*gV*(estb-ea)
    LL[LL<0]<-0
    Rabs+fluxother-Rem-LL-HH
  }
  cp <- cpair(tair)
  sb <- 5.67*10^-8
  lambda <- 42.575 * tair + 44994
  if (method != "iter") {
    gr <- (4 * em * sb * (tair + 273.15)^3) / cp
    gHr <- gH + gr
    Rem <- em * sb * (tair + 273.15)^4
    es <- satvap(tair, ice = T)
    delta<-4098*(0.6108*exp(17.27*tair/(tair+237.3)))/(tair+237.3)^2
    DD <- es-ea
    s <- delta / pk
    if (method == "open") {
      tp <- Rabs + fluxother - Rem - (surfwet*lambda*gV)*(DD/pk)
      btm <- cp * gHr + surfwet*lambda*s*gV
      OT <- tair +  tp/btm
    } else {
      lft <- (Rabs + fluxother - Rem) / (gHr * cp)
      if (surfwet > 0) {
        Gamma <- (cp * pk) / (0.622 * lambda * surfwet)
        gstar <- Gamma * (gHr/gV)
        m <- gstar / (s + gstar)
        rgt <- (es - ea) / (pk * gstar)
        OT <- tair + m * (lft - rgt)
      } else {
        OT <- tair + lft
      }
    }
  } else {
    B<-0
    x<-c(1:20)
    for (i in x) {
      Tb<-i*10-50
      B[i] <- Balance(Tb,sb,em,cp,gH,gV,lambda,surfwet,ea,pk,fluxother)
    }
    sel<-max(which(B>0))
    Tb10<-x[sel]*10-50
    B<-0
    x<-c(-2:9)
    for (i in x) {
      Tb <- Tb10+i
      B[i+3] <- Balance(Tb,sb,em,cp,gH,gV,lambda,surfwet,ea,pk,fluxother)
    }
    sel<-max(which(B>0))
    Tb1<-Tb10+x[sel]
    B<-0
    for (i in x) {
      Tb <- Tb1+i/10
      B[i+3] <- Balance(Tb,sb,em,cp,gH,gV,lambda,surfwet,ea,pk,fluxother)
    }
    sel<-max(which(B>0))
    Tb2<-Tb1+x[sel]/10
    B<-0
    for (i in x) {
      Tb <- Tb2+i/100
      B[i+3] <- Balance(Tb,sb,em,cp,gH,gV,lambda,surfwet,ea,pk,fluxother)
    }
    sel<-max(which(B>0))
    Tb3<-Tb2+x[sel]/100
    B<-0
    for (i in x) {
      Tb <- Tb3+i/1000
      B[i+3] <- Balance(Tb,sb,em,cp,gH,gV,lambda,surfwet,ea,pk,fluxother)
    }
    sel<-which(abs(B)==min(abs(B)))
    OT<-Tb3+x[sel]/1000
  }
  OT
}
