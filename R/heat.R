#' Calculate molar density of air
#' @param tc temperature (deg C)
#' @param pk atmospheric pressure (kPa)
#' @return molar density of air (mol / m3)
#' @export
phair <- function(tc = 15, pk = 101.3) {
  tk <- tc + 273.15
  ph <- 44.6 * (pk / 101.3) * (273.15 / tk)
  ph
}
#' Calculate molar specific heat of air
#' @param tc temperature (deg C)
#' @return specific heat of air at constant pressure (J / mol / K)
#' @export
cpair <- function(tc) {
  cp <-  2e-05 * tc^2 + 0.0002 * tc + 29.119
  cp
}
#' Calculate zero plane displacement
#'
#' @param hgt canopy height (m)
#' @param PAI Plant area index
#' @param cdl Control parameter scaling d/h (from Raupach 1994)
#' @return zero plane displacement (m)
#' @references Raupach, M.R. (1994) Simplified expressions for vegetation roughness length and zero-plane
#' displacement as functions of canopy height and area index. Boundary-Layer Meteorology 71: 211-216.
#' @export
zeroplanedis<-function(hgt, PAI, cdl=7.5) {
  d<-(1-(1-exp(-sqrt(cdl*PAI)))/sqrt(cdl*PAI))*hgt
  d
}
#' Calculate roughness length governing momentum transfer
#'
#' @param hgt canopy height (m)
#' @param PAI Plant area index
#' @param zm0 substrate-surface drag coefficient
#' @param cdl Control parameter scaling d/h (from Raupach 1994)
#' @param CR roughness-element drag coefficient
#' @param psih parameter characterizes the roughness sublayer depth
#' @param umx Maximum ratio of wind velocity to friction velocity
#' @return momentum roughness length (m)
#' @references Raupach, M.R. (1994) Simplified expressions for vegetation roughness length and zero-plane
#' displacement as functions of canopy height and area index. Boundary-Layer Meteorology 71: 211-216.
#' @export
roughlength<-function(hgt, PAI, zm0 = 0.003, cdl = 7.5, CR = 0.3,
                      psih = 0.193, umx = 0.3) {
  d<-zeroplanedis(hgt,PAI,cdl)
  ur<-sqrt(zm0+(CR*PAI)/2)
  ur[ur>umx]<-umx
  zm<-(hgt-d)*exp(-0.4*ur-psih)
  zm <- ifelse(zm < zm0, zm0, zm)
  zm
}
#' Calculates forced or free laminer conductance
#'
#' @param d chacteristic dimension of surface (m)
#' @param u wind speed (m/s)
#' @param tc temperature (deg C)
#' @param dtc estimate of temperature differences of surface and air, e.g.
#' from previous time step (see details)
#' @param pk atmospheric pressure (KPa). used for calculating molar density
#' of air
#' @param dtmin minimum tmeperature difference for calculating minimum conductance (see differenc)
#' @return conductance (mol / m^2 / sec)
#' @export
#'
#' @details Calculates conductance under forced and free convection and selects
#' whichever is greater. For conductance under free convection an estimate of
#' the temperature difference between the air and  surface is needed (usually
#' from the previous timestep). Because the temperature difference in the previous
#' timestep is used when conductance is under free convection, the model can de-stabalise
#' as dtc approaches 0, leading to very low conductance and very high temperature differences in the current timestep. The parameter
#' dtmin sets a minimum temperature difference, and prevents conductance being too low.
#'
#' @examples
#' # As function of `d`
#' d <- c(1:1000) / 100
#' g <- gforcedfree(d, 2, 11, 2)
#' plot(g~d, type = "l")
#' # As function of `u`
#' u <- c(1:200)/100
#' g<-gforcedfree(0.1,u,11,4)
#' plot(g~u, type = "l")
gforcedfree <- function(d, u, tc, dtc, pk = 101.3, dtmin = 1) {
  # sets temperature difference for minimum conductance
  dtc <-ifelse(dtc < dtmin, dtmin, dtc)
  # Reynolds etc
  tk <- tc + 273.15
  v <- (0.0908 * tk - 11.531) / 1000000
  Dh <- (0.1285 * tk - 16.247) / 1000000
  Gr <- (9.807  * d^3  * abs(dtc - 0)) / (tk * v^2)
  re <-  (u * d) / v
  Pr <- v / Dh
  # Forced
  ph <- phair(tc, pk)
  cp <-  cpair(tc)
  gh <- (0.34 * Dh * ph * re^0.5 * Pr ^ (1/3)) / d
  # Free
  m <- ifelse(dtc > 0, 0.54, 0.26)
  ghfr <- (m * ph * Dh * (Gr * Pr) ^ (1/4)) / d
  # Set to whichever is higher
  gh <- ifelse(gh > ghfr, gh, ghfr)
  gh
}
#' Calculate mixing length for canopy air transport
#'
#' @param hgt canopy height (m)
#' @param PAI plant area index
#' @return mixing length (m). See details
#' @export
#' @details The mixing length is us used to calculate turbulent air transport
#' inside vegetated canopies. It is made equivalent to the above canopy value
#' at the canopy surface.
mixinglength <- function(hgt, PAI) {
  d<-.zeroplanedis(hgt,PAI)
  zm<-.roughlength(hgt,PAI)
  l_m<-(0.32*(hgt-d))/log((hgt-d)/zm)
  l_m
}
#' Calculate mixing length for canopy air transport
#'
#' @description Calculates mixing length for canopy air transport - equivelent to mean
#' leaf spacing.
#'
#' @param hgt canopy height (m)
#' @param PAI plant area index
#' @param cd drag coefficient
#' @param iw relative turbulence intensity
#' @param phi_m diabatic correction factor
#' @return dimensionless attenuation coefficient for calculating the wind height profile
#' in vegetation canopies.
#' @export
attencoef <- function(hgt, PAI = 3, cd = 0.2, iw = 0.5, phi_m  = 1) {
  l_m <- mixinglength(hgt, PAI)
  a <- ((cd * PAI * hgt) / (2 * l_m * iw))^0.5
  a <- a * phi_m^0.5
  a
}
#' Calculate conductance under turbulent convection above canopy
#'
#' @param u wind speed at height `zu` (m/s)
#' @param zu height of wind speed measurement (m)
#' @param z1 upper height to which conductance is wanted (m)
#' @param z0 lower height from which conductance is wanted (set to heat
#' exchange surface of canopy if NA)
#' @param hgt height of canopy (m)
#' @param PAI plant area index for determining canopy roughness length
#' @param tc temperature used for calculating molar density of air (deg C)
#' @param psi_m diabatic correction factor for momentum transfer
#' @param psi_h diabataic correction factor for heat transfer
#' @param zm0 minimum surface roughness
#' @param pk atmospheric pressure (kPA) used for calculating molar density of air
#' @return molar conductance under turbulent convection above canopy (mol/m^2/s)
#' @export
gturb <- function(u, zu = 2, z1, z0 = NA, hgt, PAI= 3, tc = 15,
                  psi_m = 0, psi_h = 0, zm0 = 0.004, pk = 101.3) {
  d <- zeroplanedis(hgt, PAI)
  zm <- roughlength(hgt, PAI, zm0)
  zh <- 0.2 * zm
  if (is.na(z0)) z0 <- d + zh
  ph <- phair(tc, pk)
  cp <- cpair(tc)
  lnr<-log((z1-d)/zh)/log((zu-d)/zh)-log((z0-d)/zh)/log((zu-d)/zh)
  ustr <- (0.4 * u) / (log((zu - d) / zm) + psi_m)
  g <- (0.4 * ph * ustr) / (log((z1 - d) / (z0 - d)) + psi_h*lnr)
  Dh <- (0.1285 * (tc + 273.15) - 16.247) / 1000000
  zref <- d + zh
  gmin <- (ph * Dh) / (z1 - z0)
  g <- ifelse(g < gmin, gmin, g)
  g[is.na(g)] <- gmin
  g
}
#' Calculates conductance under turbulent convection within the canopy
#'
#' @param uh wind speed at height of canopy top (m/s) as returned by [microclimc::windprofile()]
#' @param z1 upper height to which conductance is wanted (m)
#' @param z0 lower height from which conductance is wanted (m)
#' @param tc1 temperature of upper layer (dec C) usually in previous timestep
#' @param tc0 temperature of lower layer (dec C) usually in previous timestep
#' @param hgt height of canopy (m)
#' @param PAI plant area index for mixing and attenuation coefficients
#' @param cd drag coefficient
#' @param iw relative turbulence intensity
#' @param phi_m diabatic correction factor
#' @param pk atmospheric pressure used for calculating moar density of air (kPA)
#' @return conductance under tubulent convection within the canopy (mol/m^2/sec)
#' @export
gcanopy <- function(uh, z1, z0, tc1, tc0, hgt, PAI = 3, cd = 0.2,
                    iw = 0.5, phi_m = 1, pk = 101.3) {
  a <- attencoef(hgt, PAI, cd, iw, phi_m)
  l_m <- mixinglength(hgt, PAI)
  tcm <- (tc1 + tc0) / 2
  ph <- phair(tcm, pk)
  tp <- (l_m * iw * ph * uh * a) / phi_m
  e0 <- exp(-a * (z0 / hgt - 1))
  e1 <- exp(-a * (z1 / hgt - 1))
  g <- tp / (e0 - e1)
  # Set minimum
  cp <- cpair(tcm)
  dTdz <- abs((tc1 - tc0) / (z1 - z0))
  Kmin <- (iw * l_m^2 * 3 / 0.74)^ (2/3) *  ((4.6 * dTdz) / (tcm + 273.15))^0.5
  gmin <- (Kmin * ph) / (z1 - z0)
  g <- ifelse(g < gmin, gmin, g)
  g[is.na(g)] <- mean(gmin,na.rm=T)
  g
}
#' Calculates diabatic correction factor above canopy
#'
#' @description Calculates the diabatic correction factors used in adjustment of
#' wind profiles and calculation of turbulent conductivity above canopy
#'
#' @param tc temperature
#' @param pk atmospheric pressure (kPa)
#' @param H Heat flux (W / m^2)
#' @param uf friction velocity (m/s)
#' @param zi height to which correction factor is wanted (m)
#' @param d zero plane displacement height as returned by [zeroplanedis()]
#' @return a list with the following components:
#' @return `psi_m` diabatic correction factor for momentum transfer
#' @return `psi_h` diabatic correction factor for heat transfer
#' @seealso [diabatic.approx()]
#' @export
diabatic_cor <- function(tc, pk = 101.3, H = 0, uf, zi = 2, d) {
  Tk <- tc + 273.15
  ph <- phair(tc, pk)
  cp <-  cpair(tc)
  st <- -(0.4 * 9.81 * (zi - d) * H) / (ph * cp * Tk * uf^3)
  # Stable flow
  sel <- which(st < 0) # unstable
  # Stable
  psi_h <- suppressWarnings(6 * log(1 + st))
  psi_m <- psi_h
  # Unstable
  psi_h[sel] <-   -2 * log((1 + (1 - 16 * st[sel])^0.5) / 2)
  psi_m[sel] <- 0.6 * psi_h[sel]
  psi_m <-ifelse(psi_m > 5, 5, psi_m)
  psi_h <-ifelse(psi_h > 5, 5, psi_h)
  return(list(psi_m = psi_m, psi_h = psi_h))
}
#' Calculates approximation of diabatic correction factor above canopy
#'
#' @description Calculates approximation of diabatic correction factor above canopy without
#' using values in current time-step only.
#' @param tc temperature (deg C)
#' @param u wind speed at 2 m height
#' @param H Heat flux (W / m^2)
#' @param hgt hgt of canopy (m)
#' @param zi height to which correction factor is wanted (m)
#' @return a list with the following components:
#' @return `psi_m` diabatic correction factor for momentum transfer
#' @return `psi_h` diabatic correction factor for heat transfer
#' @details The function `diabatic.approx` bi-passes the need to run microclimate
#' models in hourly time-steps
#' @seealso [diabatic_cor()] [diabatic_cor_can()]
#' @export
diabatic.approx <- function(tc, u, H, hgt, PAI, zi) {
  diabfu<-function(tc,u,H,hgt,PAI,zi) {
    wfa<-(0.0527*u^2-1.6741*u-5.4916)/(-7.113)
    wfa[wfa<1]<-1
    lwfb<-(0.1292*u+0.3979)/0.5271
    a<-(1.9932*log(hgt)-8.046)*wfa
    lb<- (-0.056*log(hgt)+0.2614)*lwfb
    b<-1/(1+exp(-lb))
    d<-zeroplanedis(hgt, PAI)
    zm<-roughlength(hgt, PAI)
    uf<-(u/0.4)/log(2/zm)
    Hf<-H/(uf^3*(tc+273.15))
    pshu<- -2*log((1+(1-a*Hf^b)^0.5)/2)
    psmu<- 0.6*pshu
    return(list(psi_m=psmu,psi_h=pshu))
  }
  diabfs<-function(tc,u,H,hgt,PAI,zi) {
    p1<-6*log(1+3.733903*Hf^1.247703)
    p2<-6*log(1+25.61992*Hf^1.88255)
    p3<-6*log(1+7195*Hf^4.04)
    wgt1<-1/abs(p1-0.25)
    wgt2<-1/abs(p1-0.75)
    wgt3<-1/abs(p1-1.5)
    ws1<-wgt1+wgt2
    wgt1<-wgt1/ws1
    wgt2<-wgt2/ws1
    pred1<-ifelse(p1<0.25,p1,wgt1*p1+wgt2*p2)
    ws2<-wgt2+wgt3
    wgt2<-wgt2/ws2
    wgt3<-wgt3/ws2
    pred2<-ifelse(p1<=1,wgt2*p2+wgt3*p3,p3)
    pred<-ifelse(p1<=0.75,pred1,pred2)
    pred<-pred*0.927
    Pmult<- -0.043*log(PAI)+1.0018
    Hmult<- -0.6168*hgt+1.2995
    pshs<-pred*Pmult*Hmult
    pshs<-ifelse(pshs>3,3,pshs)
    psms<-pshs
    return(list(psi_m=psms,psi_h=pshs))
  }
  diabf2<-function(tc,u,H,hgt,PAI,zi) {
    db1<-diabfu(tc,u,H,hgt,PAI,zi)
    db2<-diabfs(tc,u,H,hgt,PAI,zi)
    psi_m<-db1$psi_m
    psi_h<-db1$psi_h
    sel<-which(H<0)
    psi_m[sel]<-db2$psi_m[sel]
    psi_h[sel]<-db2$psi_h[sel]
    rath<-log((zi-d)/(0.2*zm))/log((2-d)/(0.2*zm))
    ratm<-log((zi-d)/zm)/log((2-d)/zm)
    psi_m<-psim*ratm
    psi_h<-psih*rath
    return(list(psi_m=psi_m,psi_h=psi_h))
  }
  db<-suppressWarnings(diabf2(tc,u,H,hgt,PAI,zi))
  psi_m<-db$psi_m
  psi_h<-db$psi_h
  sel<-which(hgt>2)
  psi_m[sel]<-0
  psi_h[sel]<-0
  return(list(psi_m=psi_m,psi_h=psi_h))
}
#' Calculates diabatic correction factor in canopy
#'
#' @description Calculates the diabatic correction factors used in adjustment of
#' wind profiles and calculation of turbulent conductivity within the canopy
#'
#' @param tc vector of air temperatures at canopy nodes (deg C)
#' @param uz vector of wind speeds at canopy nodes (m/s)
#' @param z vector of heights of canopy nodes (m)
#' @param PAI vector of plant area index values at each canopy node
#' @return a list with the following components:
#' @return `phi_m` diabatic correction factor for momentum transfer
#' @return `phi_h` diabatic correction factor for heat transfer
#'
#' @export
diabatic_cor_can <- function(tc, uz, z, PAI) {
  dtc<-tc[2:length(tc)]-tc[1:(length(tc)-1)]
  dz<-z[2:length(z)]-z[1:(length(z)-1)]
  dtdz<-dtc/dz
  tk <- (tc[2:length(tc)]+tc[1:(length(tc)-1)])/2+273.15
  PAIl<-(PAI[2:length(PAI)]+PAI[1:(length(PAI)-1)])/2
  l_m <- mixinglength(dz, PAIl)
  u<-(uz[2:length(uz)]+uz[1:(length(uz)-1)])/2
  Ri <- (9.81 / tk) * dtdz * (l_m / u)^2
  Ri[Ri > 0.15] <- 0.15
  Ri[Ri <= -0.1120323] <- -0.112032
  st <- (0.74 * (1 + 8.926 * Ri) ^ 0.5 + 2 * 4.7 * Ri - 0.74) / (2 * 4.7 * (1 - 4.7 * Ri))
  sel <- which(dtdz <= 0) # unstable
  st[sel]<-Ri[sel]
  # Stable
  phi_m <- 1 + (6 * st) / (1 + st)
  phi_h <- phi_m
  # Unstable
  phi_m[sel] <- 1 / (1 - 16 * st[sel])^0.25
  phi_h[sel] <- phi_m[sel]^2
  return(list(phi_m = mean(phi_m), phi_h = mean(phi_h)))
}
