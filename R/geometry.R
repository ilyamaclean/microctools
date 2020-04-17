globalVariables(c("globalclimate", "habitats", "soilparams", "climvars", "weather"))
#' Generates plant area index profile
#'
#' @param m number of canopy nodes
#' @param PAI total plant area index for canopy
#' @param skew number between 0 and 10 indicating the degree of skew towards top of
#' canopy in canopy foliage (see details)
#' @param spread positive non-zero number less than 100 indicating the degree of spread in
#' canopy foliage (see details)
#' @details in specifyign `skew`, lower numbers indicate greater skew towards top of
#' canopy (5 = symetrical). In specifying `spread` a value of one indicates almost
#' all the foliage in concentrated in one canopy layer, whereas a value of 100 indicates
#' completely evenly spread.
#' @return a vector of length `m` of plant area indices for each canopy layer
#' @export
#' @importFrom stats dbeta
#' @examples
#' pai <- PAIgeometry(100, 10, 7, 70)
#' plot(pai, type = "l")
PAIgeometry <- function(m, PAI, skew, spread) {
  skew <- 10 - skew
  # Plant area index of canopy layer
  shape1 <- 100 / spread
  x<-c(1:m)/(m+1)
  if (skew > 5) {
    shape2 <- (10 - skew) / 5 + 1
    shape2 <- shape2 / 2 * shape1
    y <- rev(dbeta(x, shape1, shape2))
  } else {
    shape2 <- (skew + 5) / 5
    shape2 <- shape2 / 2 * shape1
    y <- dbeta(x, shape1, shape2)
  }
  y <- PAI / sum(y) * y
  y
}
#' Generates vegetation thickness profile
#'
#' @description `thickgeometry` is used is to calculate the mean thickness of leaves
#' and woody vegetation such that volume = Plant Area Index * thickness
#'
#' @param m number of canopy nodes
#' @param mthick mean thickness of vegetation
#' @param thmax mean thickness of lowest portion of canopy
#' @param thmin mean thickness of upper portion of canopy
#' @return a vector of thicknesses for each canopy node (m)
#' @export
#' @examples
#' pai <- PAIgeometry(1000, 10, 7, 70)
#' thick <- thickgeometry(1000,0.4,0.7,0.1)
#' hgt <- 25
#' dens <- thick * pai * length(pai) / hgt
#' z <- c(1:1000) / 40
#' par(mar=c(5, 5, 2, 2))
#' plot(z ~ dens, type = "l",  cex.axis = 1.5, cex.lab = 1.5,
#'     xlab = expression(paste("Volume ",(~m^3 / m^3))),
#'     ylab = "Height (m)")
thickgeometry <- function(m, mthick, thmax, thmin) {
  x <- c(1:m)/m
  me <- 0
  for (i in 1:1000) {
    y<- (thmax - thmin) * x^(i/80) + thmin
    me[i] <- mean(y)
  }
  dif <- abs(me-mthick)
  sel <- which(dif == min(dif))
  y<- (thmax - thmin) * x^(sel/80) + thmin
  y <- rev(y)
  y

}
#' Generates fraction of green leaves profile
#'
#' @param m number of canopy nodes
#' @param pLAI mean proportion of green leaf foliage
#' @param skew number between 0 and 10 indicating the degree of skew towards top of
#' canopy in green foliage (see details).
#' @return vector of length `m` of the proportion of plant area that is green leaf area
#' @export
#' @examples
#' plot(LAIfrac(100,0.8, 0),type="l", ylim = c(0,1))
#' plot(LAIfrac(100,0.8, 5),type="l", ylim = c(0,1))
#' plot(LAIfrac(100,0.8, 10),type="l", ylim = c(0,1))
LAIfrac <- function(m, pLAI, skew = 7) {
  skew <- skew/2 + 5
  x <- c(1:m)/m
  a <- 10^7/(10^skew)
  rge <- (1/a)
  mL <- 0
  for (i in 1:1000) {
    mu <- (i/100) * rge
    y1 <- a + x
    y1 <- mu * y1
    y <- 2 / (1+exp(-y1)) - 1
    mL[i] <- mean(y)
  }
  dif <- abs(mL - pLAI)
  sel <- which(dif == min(dif))
  mu <- (sel/100) * rge
  y1 <- a + 2*x
  y1 <- mu * y1
  ipLAI <- 2 / (1+exp(-y1)) - 1
  ipLAI
}
#' Generates relative turbulence intensity profile
#'
#' @param m number of canopy nodes
#' @param iwmin minimum relative turbulence intensity
#' @param iwmax maximum relative turbulence intensity
#' @param increasing logical indicating whether turbulence intensity increases (TRUE) or decreases
#' (FALSE) with height
#' @return relative turbulence intensity for each node `m`.
#' @details default values are for a maize crop Shaw et al (1974) Agricultural Meteorology,
#' 13: 419-425.
#' @export
iwgeometry <- function(m, iwmin = 0.36, iwmax = 0.9, increasing = FALSE) {
  x <- c(1:m) / m
  if (increasing) {
    y <- iwmin + (iwmax - iwmin) * x
  } else y <- iwmax - (iwmax - iwmin) * x
  y
}
#' Derives leaf area index, leaf geometry and canopy height from habitat
#'
#' @description `PAIfromhabitat` generates an hourly dataset for an entire year of
#' leaf area index values, the ratio of vertical to horizontal projections of leaf
#' foliage and canopy height from habitat.
#'
#' @param habitat a character string or numeric value indicating the habitat type. See [habitats()]
#' @param lat a single numeric value representing the mean latitude of the location for which the solar index is required (decimal degrees, -ve south of the equator).
#' @param long a single numeric value representing the mean longitude of the location for which the solar index is required (decimal degrees, -ve west of Greenwich meridian).
#' @param meantemp an optional numeric value of mean annual temperature (ºC) at reference height.
#' @param cvtemp an optional numeric value of the coefficient of variation in temperature (K per 0.25 days) at reference height.
#' @param rainfall an optional numeric value mean annual rainfall (mm per year)
#' @param cvrain an optional numeric value of the coefficient of variation in raifall (mm per 0.25 days) at reference height.
#' @param wetmonth an optional numeric value indicating which month is wettest (1-12).
#' @param year the year for which data are required.

#' @return a list with the following items:
#' @return `lai` hourly leaf area index values,
#' @return `x` the ratio of vertical to horizontal projections of leaf foliage
#' @return `height` the heigbht of the canopy in metres
#' @return `obs_time` an object of class POSIXlt of dates and times coressponding to
#' each value of `lai`
#' @import sp raster
#' @importFrom stats spline
#' @export
#'
#' @details
#' If no values of `meantemp`, `cvtemp`, `rainfall`, `cvrain` or `wetmonth` are
#' provided, values are obtained from [globalclimate()]. Variable `lai` is derived by fitting
#' a Gaussian curve parameters of which are climate and location-dependent. Functional
#' fits were calibrated using MODIS data (see https://modis.gsfc.nasa.gov/data/).
#'
#' @examples
#' pxh <- PAIfromhabitat("Deciduous broadleaf forest", 50, -5.2, 2015)
#' pxh$height
#' pxh$x
#' plot(pxh$lai ~ as.POSIXct(pxh$obs_time), type = "l", xlab = "Month",
#'      ylab = "LAI", ylim = c(0, 6))

PAIfromhabitat <- function(habitat, lat, long, year, meantemp = NA, cvtemp = NA,
                           rainfall = NA, cvrain = NA, wetmonth = NA) {
  laigaus <- function(minlai, maxlai, pkday, dhalf, yr) {
    diy <- 365
    sdev <- 0.0082 * dhalf^2 + 0.0717 * dhalf + 13.285
    difv <- maxlai - minlai
    x<-c(-diy:diy)
    y <- 1 / (sdev * sqrt(2 * pi)) * exp(-0.5 * (((x - 0) / sdev) ^ 2))
    y[(diy + ceiling(0.5 * diy)):(2 * diy + 1)] <- y[(diy - ceiling(0.5 * diy)):diy]
    st <- diy + 1 - pkday
    y <- y[st:(st + diy - 1)]
    x <- c(1:diy)
    x <- c(0, x, c(366:375))
    y <- c(y[diy], y, y[1:10])
    sel <-c(0:15) * 25 + 1
    x<-x[sel]
    y<-y[sel]
    tme <- as.POSIXct((x * 24 * 3600), origin = paste0(yr - 1,"-12-31 12:00"), tz = "GMT")
    xy <- spline(tme, y, n = diy * 24 + 241)
    tme2 <- as.POSIXlt(xy$x, origin = "1970-01-01 00:00", tz = "GMT")
    sel <- which(tme2$year + 1900 == yr)
    y <- xy$y[sel]
    dify <- max(y) - min(y)
    y <- y * (difv / dify)
    y <- y + minlai - min(y)
    return(y)
  }
  long <- ifelse(long > 180.9375, long - 360, long)
  long <- ifelse(long < -179.0625, long + 360, long)
  ll <- SpatialPoints(data.frame(x = long, y = lat))
  diy <- 366
  if (year%%4 == 0) diy <- 366
  if (year%%100 == 0 & year%%400 != 0) diy <- 365
  mmonth <-c(16, 45.5, 75, 105.5, 136, 166.5, 197, 228, 258.5, 289, 319.5, 350)
  if (diy == 365) mmonth[2:12] <- mmonth[2:12] + 0.5
  e <- extent(c(-179.0625, 180.9375, -89.49406, 89.49406))
  clim <- c(meantemp, cvtemp, rainfall, cvrain, wetmonth)
  for (i in 1:5) {
    if (is.na(clim[i])) {
      r <- raster(globalclimate[,,i])
      extent(r) <- e
      clim[i] <- extract(r, ll)
    }
  }
  wgts <- function(x1, x2, ll, lmn, lmx) {
    ll <- ifelse(ll < lmn, lmn, lat)
    ll <- ifelse(ll > lmx, lmx, lat)
    w <- 1 - (abs(ll - lmn)  / (abs(ll - lmn)  + abs(ll - lmx)))
    y <- w * x1 + (1 - w) * x2
    y
  }
  # By habitat type
  if (habitat == "Evergreen needleleaf forest" | habitat == 1) {
    h2 <- 74.02 + 5.35 * clim[1]
    h1 <-  203.22 - 35.63 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 50, 50, hperiod)
    p2 <- 216.71 - 2.65 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <- 2.33 + 0.0132 * clim[1]
    } else maxlai <- 2.33 + 0.0132 * 20
    minlai <- 1.01
    x <- 0.4
    hgt <- 15
  }
  if (habitat == "Evergreen Broadleaf forest" | habitat == 2) {
    hperiod <-  154.505 + 2.040 * clim[1]
    hperiod <- ifelse(hperiod < 50, 50, hperiod)
    peakdoy <- peakdoy <- mmonth[round(clim[5], 0)]
    maxlai <- 1.83 + 0.22 * log(clim[3])
    minlai <- (-1.09) + 0.4030 * log(clim[3])
    minlai <- ifelse(minlai < 1, 1, minlai)
    x <- 1.2
    hgt <- 20
  }
  if (habitat == "Deciduous needleleaf forest" | habitat == 3) {
    h2 <- 51.18 + 3.77  * clim[1]
    h1 <- 152
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    p2 <- 204.97 - 1.08 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <-  2.62 + 0.05 * clim[1]
    } else maxlai <- 2.62 + 0.05 * 20
    minlai <- 0.39
    x <- 0.4
    hgt <- 10
  }
  if (habitat == "Deciduous broadleaf forest" | habitat == 4) {
    h2 <- 47.6380 + 2.9232 * clim[1]
    h1 <- 220.06 - 79.19 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 32.5, 32.5, hperiod)
    p2 <- 209.760 - 1.208 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <- 3.98795 + 0.03330 * clim[1]
    } else maxlai <- 3.98795 * 0.03330 * 20
    minlai <- 0.4808
    x <- 1.2
    hgt <- 15
  }
  if (habitat == "Mixed forest" | habitat == 5) {
    h2 <- 74.02 + 5.35 * clim[1]
    h1 <-  203.22 - 35.63 * clim[4]
    hperiod1 <- wgts(h1, h2, abs(lat), 0, 20)
    h2 <- 51.18 +  3.77  * clim[1]
    h1 <-  152
    hperiod2 <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- (hperiod1 + hperiod2) / 2
    hperiod <- ifelse(hperiod < 30.5, 30.5, hperiod)
    p2 <- 216.71 - 2.65 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy1 <- wgts(p1, p2, abs(lat), 0, 30)
    p2 <-  204.97 - -1.08 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    peakdoy2 <- wgts(p1, p2, abs(lat), 0, 30)
    peakdoy <- (peakdoy1 + peakdoy2) / 2
    if (clim[1] <= 20) {
      maxlai1 <- 2.33 + 0.0132 * clim[1]
      maxlai2 <-  2.62 + 0.05 * clim[1]
      maxlai <- (maxlai1 + maxlai2) / 2
    } else maxlai <- 3.107
    minlai <- 0.7
    x <- 0.8
    hgt <- 10
  }
  if (habitat == "Closed shrublands" | habitat == 6) {
    h2 <- 33.867 + 6.324 * clim[1]
    h1 <-  284.20 - 102.51 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 30.5, 30.5, hperiod)
    p2 <- 223.55 - 3.125 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    maxlai <- 2.34
    minlai <- -0.4790 + 0.1450 * log(clim[3])
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
    x <- 1
    hgt <- 2
  }
  if (habitat == "Open shrublands" | habitat == 7) {
    h2 <- 8.908 + 4.907 * clim[1]
    h1 <-  210.09 - 28.62 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 38.3, 38.3, hperiod)
    p2 <- 211.7 - 4.085 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    maxlai <- -0.7206 + 0.272 * log(clim[3])
    minlai <- -0.146 +  0.059 * log(clim[3])
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
    x <- 0.7
    hgt <- 1.5
  }
  if (habitat == "Woody savannas" | habitat == 8) {
    hperiod1 <-  47.6380 + 2.9232 * clim[1]
    hperiod1 <- ifelse(hperiod1 < 32.5, 32.5, hperiod1)
    hperiod2 <- 71.72 + 3.012 * clim[1]
    h1 <- (hperiod1 + hperiod2) / 2
    h2 <- 282.04 - 92.28 * clim[4]
    h2 <- ifelse(hperiod1 < 31.9, 31.9, hperiod1)
    hperiod <- wgts(h1, h2, abs(lat), 25, 35)
    peakdoy1 <- 209.760 - 1.208 * clim[1]
    peakdoy1 <- ifelse(peakdoy1 > 244, 244, peakdoy1)
    if (lat < 0)  peakdoy1 <- ( peakdoy1 + diy / 2)%%diy
    peakdoy2 <- 211.98 - 3.4371 * clim[1]
    peakdoy2 <- ifelse(peakdoy2 > 244, 244, peakdoy2)
    if (lat < 0)  peakdoy2 <- (peakdoy2 + diy / 2)%%diy
    p2 <- (peakdoy1 + peakdoy2) / 2
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 40)
    if (clim[1] <= 20) {
      maxlai1 <- 3.98795 + 0.03330 * clim[1]
      maxlai2 <- 1.0532 * 0.016 * clim[1]
    } else {
      maxlai1 <- 3.98795 + 0.03330 * 20
      maxlai2 <- 1.0532 * 0.016 * 20
    }
    mx2 <- (maxlai1 + maxlai2) / 2
    minlai1 <- 0.4808
    minlai2 <- 0.0725 * 0.011 * clim[1]
    mn2 <- (minlai1 + minlai2) / 2
    mx1 <- 1.298 + 0.171 * log(clim[3])
    mn1 <- -2.9458 + 0.5889 * log(clim[3])
    maxlai <- wgts(mx1, mx2, abs(lat), 10, 40)
    minlai <- wgts(mn1, mn2, abs(lat), 10, 40)
    minlai <- ifelse(minlai < 0.0362, 0.0362, minlai)
    x <- 0.7
    hgt <- 3
  }
  if (habitat == "Savannas" | habitat == 9 |
      habitat == "Short grasslands" | habitat == 10 |
      habitat == "Tall grasslands" | habitat == 11) {
    h2 <- 71.72 + 3.012 * clim[1]
    h1 <- 269.22 -  89.79 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 31.9, 31.9, hperiod)
    p2 <- 211.98 - 3.4371 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    mx2 <- 1.48
    mn2 <- 0.0725 * 0.011 * clim[1]
    mx1 <- 0.1215 + 0.2662 * log(clim[3])
    mn1 <- 0.331 + 0.0575 * log(clim[3])
    maxlai <- wgts(mx1, mx2, abs(lat), 10, 40)
    minlai <- wgts(mn1, mn2, abs(lat), 10, 40)
    minlai <- ifelse(minlai < 0.762, 0.762, minlai)
    x <- 0.15
    if (habitat == "Savannas" | habitat == 9) hgt <- 1.5
    if (habitat == "Short grasslands" | habitat == 10) hgt <- 0.25
    if (habitat == "Tall grasslands" | habitat == 11) hgt <- 1.5
  }
  if (habitat == "Permanent wetlands" | habitat == 12) {
    h2 <- 76 + 4.617 * clim[1]
    h1 <- 246.68 - 66.82 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 40, 40, hperiod)
    p2 <- 219.64 - 2.793 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    maxlai <- -0.1782 + 0.2608 * log(clim[3])
    maxlai <- ifelse(maxlai < 1.12, 1.12, maxlai)
    minlai <-  -0.1450 + 0.1440 * log(clim[3])
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
    x <- 1.4
    hgt <- 0.5
  }
  if (habitat == "Croplands" | habitat == 13) {
    h2 <- 54.893 +  1.785 * clim[1]
    h1 <- 243 - 112.18 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 10, 30)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 212.95 - 5.627 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <- 3.124 - 0.0886 * clim[1]
    } else maxlai <- 3.124 - 0.0886 * 20
    maxlai <- ifelse(maxlai > 3.14, 3.14, maxlai)
    minlai <- 0.13
    x <- 0.2
    hgt <- 0.5
  }
  if (habitat == "Urban and built-up" | habitat == 14) {
    h2 <- 66.669 +  5.618 * clim[1]
    h1 <- 283.44 - 86.11 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 215.998 - 4.2806 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 30)
    if (clim[1] <= 20) {
      maxlai <- 1.135 - 0.0244 * clim[1]
    } else maxlai <- 1.135 - 0.0244 * 20
    maxlai <- ifelse(maxlai > 1.15, 1.15, maxlai)
    minlai <- 0.28
    x <- 1
    hgt <- 0.8
  }
  if (habitat == "Cropland/Natural vegetation mosaic" | habitat == 15) {
    h2 <- 29.490 +  8.260 * clim[1]
    h1 <- 326.46 - 161.70 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 10, 30)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 210.867 - 3.5464 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 30)
    if (clim[1] <= 20) {
      maxlai <- 3.5485 - 0.09481 * clim[1]
    } else maxlai <- 3.5485 - 0.09481 * 20
    maxlai <- ifelse(maxlai > 3.14, 3.14, maxlai)
    if (clim[1] <= 20) {
      minlai <- -0.072815 - 0.044546 * clim[1]
    } else minlai <- -0.072815 - 0.044546 * 20
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
    x <- 0.5
    hgt <- 1
  }
  if (habitat == "Barren or sparsely vegetated" | habitat == 16) {
    h2 <- 80.557 +  6.440 * clim[1]
    h1 <- 344.65 -  -191.94 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 236.0143 - 3.4726 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 30)
    maxlai <- -0.05491 + 0.05991 * log(clim[4])
    maxlai <- ifelse(maxlai < 0.81, 0.81, maxlai)
    minlai <- 0.08
    x <- 0.6
    hgt <- 0.15
  }
  if (habitat == "Open water" | habitat == 17) {
    hperiod <- 100
    peakdoy <- 50
    maxlai <- 0
    minlai <- 0
    hgt <- 0
  }
  lai <- laigaus(minlai, maxlai, peakdoy, hperiod, year)
  if (habitat ==  "Short grasslands" | habitat == 10) lai <- lai / 2
  tme <- c(1:length(lai)) - 1
  tme <- as.POSIXlt(tme * 3600, origin = paste0(year,"-01-01 00:00"), tz = "UTC")
  return(list(lai = lai, x = x, height = hgt, obs_time = tme))
}
#' Converts daily times to hourly
#' @export
.tme.sort <- function(tme) {
  tme<-as.numeric(tme)
  tme<-as.POSIXlt(tme, origin = "1970-01-01 00:00", tz = "UTC")
  zero <- function(x) ifelse(x < 10, paste0("0", x), paste0("", x))
  dstart<-paste0(tme$year[1]+1900,"-",zero(tme$mon[1]+1),"-",zero(tme$mday[1])," 00:00")
  dfinish<-paste0(tme$year[length(tme)]+1900,"-",zero(tme$mon[length(tme)]+1),"-",
                  zero(tme$mday[length(tme)])," 23:00")
  tmeh<-seq(as.POSIXlt(dstart, format = "%Y-%m-%d %H:%M", origin = "1900-01-01", tz = "UTC"),
            as.POSIXlt(dfinish, format = "%Y-%m-%d %H:%M", origin = "1900-01-01", tz = "UTC"),
            by = 'hours')
  as.POSIXlt(tmeh)
}
#' Returns PAI for multiple years
#' @export
.PAI.sort <- function(habitat, lat, long, tme) {
  tme<-.tme.sort(tme)
  lai<-0; ota<-as.POSIXlt(0,origin="1970-01-01 00:00",tz="UTC")
  yrs<-unique(tme$year) +1900
  for (i in 1:length(yrs)) {
    pai <- PAIfromhabitat(habitat, lat, long, yrs[i])
    sel<-which(tme$year+1900 == yrs[i])
    tme1<-tme[sel]
    tmemn<-as.numeric(min(tme1))
    tmemx<-as.numeric(max(tme1))
    ot<-as.numeric(pai$obs_time)
    sel <- which(ot >= tmemn & ot <= tmemx)
    l <- pai$lai[sel]
    otx<- pai$obs_time[sel]
    lai<-c(lai,l)
    ota<-c(ota,otx)
  }
  lai<-lai[-1]; ota<-ota[-1]
  ota<-as.POSIXlt(as.numeric(ota),origin="1970-01-01 00:00",tz="UTC")
  pai$lai<-lai
  pai$obs_time<-ota
  pai
}
#' generates PAI and pLAI for a given habitat type
#' @export
.habgen<-function(habitat,lat,long,tme,m,m2,skew,spread,pLAIo1,under=TRUE) {
  # PAI
  pai <- .PAI.sort(habitat, lat, long, tme)
  m2 <- round((1 / pai$height) * m, 0)
  wgt <- (m2 / m) * 0.25
  if (class(PAIt) != "logical") {
    if (length(PAIt) == 1) {
      PAI<-PAIgeometry(m, PAIt * (1 - wgt), skew, spread)
      if (under) {
        PAIu <- PAIgeometry(m2, PAIt * wgt, 1, 50)
        PAI <- PAI + PAIu
      }
    }
    if (length(PAIt) == m) {
      mtme<-mean(as.numeric(tme))
      tdif<-abs(as.numeric(pai$obs_time)-mtme)
      sel<-which(tdif==min(tdif))
      rlai<-pai$lai[sel]
      mult<-pai$lai/rlai
      PAI <- array(NA, dim = c(m,length(pai$lai)))
      for (i in 1:m) PAI[i,] <- PAIt[i]*mult
      PAIs<-apply(PAI)
    }
    if (length(PAIt) != 1 & length(PAIt) != m) {
      tmeh<-.tme.sort(tme)
      x<-c(1:length(PAIt))
      PAIt<-spline(x,PAIt,n=length(tmeh))$y
      PAIo<-PAIgeometry(m, max(PAIt) * (1 - wgt), skew, spread)
      if (under) {
        PAIu <- PAIgeometry(m2, max(PAIt) * wgt, 1, 50)
        PAIo <- PAIo + PAIu
      }
      mn<-min(PAIo)
      PAIo<-PAIo-mn
      mu<-PAIt/max(PAIt)
      PAI <- array(NA, dim = c(m,length(pai$lai)))
      for (i in 1:m) {
        PAI[i,]<- PAIo[i]*mu+mn
      }
    }
  } else {
    xx<-PAIfromhabitat(1, lat, long, tme$year[length(tme)]+1900)$lai
    mxPAI <- max(xx,pai$lai)
    PAIo <- PAIgeometry(m, mxPAI * (1 - wgt), 7.5, 70)
    if (under) {
      PAIu <- PAIgeometry(m2, mxPAI * wgt, 1, 50)
      PAIu <- c(PAIu, rep(0, m - m2))
      PAIo <- PAIo + PAIu
    }
    if (length(tme) > 1) {
      mn<-min(PAIo)
      PAIo<-PAIo-mn
      mu<-pai$lai/mxPAI
      PAI <- array(NA, dim = c(m,length(pai$lai)))
      for (i in 1:m) {
        PAI[i,]<- PAIo[i]*mu+mn
      }
    } else PAI <- PAIo
  }
  # pLAI
  pLAIo <- LAIfrac(m, pLAIo1, 6)
  if (under) pLAIu <- c(LAIfrac(m2, 0.9, 6), rep(0, m - m2))
  if (class(PAI) == "matrix") {
    pLAI <- PAI*NA
    PAIs<-apply(PAI,2,sum)
    mxPAI<-max(PAIs)
    mxpLAI <- PAIs / mxPAI
    for (i in 1:m) {
      if (under) {
        pLAI2 <- pLAIo[i] * ((PAI[i,] - PAIu[i]) / PAI[i,]) + pLAIu[i] * (PAIu[i] / PAI[i,])
        pLAI2[is.na(pLAI2)] <- pLAIo[is.na(pLAI2)]
      } else pLAI2<-pLAIo
      pLAI[i,] <- (mxpLAI / max(pLAI2)) * pLAI2
    }
  } else {
    mxpLAI <- max(pai$lai) / mxPAI
    if (under) {
      pLAI2 <- pLAIo * ((PAI - PAIu) / PAI) + pLAIu * (PAIu / PAI)
      pLAI2[is.na(pLAI2)] <- pLAIo[is.na(pLAI2)]
    } else pLAI2<-pLAIo
    pLAI <- (mxpLAI / max(pLAI2)) * pLAI2
  }
  return(list(PAI=PAI,pLAI=pLAI))
}
#' Derives vegetation paramaters from habitat type
#'
#' @description `habitatvars` generates vegetation parameters required to run the model
#' from habitat type.
#'
#' @param habitat a character string or numeric value indicating the habitat type. See [habitats()]
#' @param lat a single numeric value representing the mean latitude of the location for which the solar index is required (decimal degrees, -ve south of the equator).
#' @param long a single numeric value representing the mean longitude of the location for which the solar index is required (decimal degrees, -ve west of Greenwich meridian).
#' @param tme a single value or vector of POSIXlt objects of date(s) and time(s). See details
#' @param m number of canopy nodes
#' @param PAIt optional single numeric value or vector of values of total plant area index for canopy (see details).
#' @return a list with the following components:
#' @return `hgt` height of vegetation (m)
#' @return `PAI` a vector or array of Plant Area Index values (see details)
#' @return `x` the ratio of vertical to horizontal projections of leaf foliage
#' @return `lw` mean leaf width (m)
#' @return `cd` drag coefficient
#' @return `iw` a vector of length `m` of relative turbulence intensities
#' @return `hgtg` height of understory vegetation (m)
#' @return `zm0` roughnes length governing momentum transfer of understory vegetation (m)
#' @return `pLAI` a vector or array of proportions of green vegetation (see details)
#' @return `refls` reflectivity (shortwave radiation) of leaves
#' @return `refg` reflectivity (shortwave radiation) of ground
#' @return `refw` reflectivity (shortwave radiation) of dead vegetation and woody stems
#' @return `reflp` #  Reflectivity of leaves to photosynthetically active radiation
#' @return `vegem` thermal emissivity of vegetation
#' @return `gsmax` maximum stomatal conductance (mol / m^2 / s)
#' @return `q50` Value of PAR when photosynthetically active radiation is at 50% of
#'          maximum (micromoles / m^2 / s^1)
#' @return `thickw` a vector of length `m` of mean thickness (m) of vegetation such that `thickw`
#'          x PAI = volume of vegetation (m^3 / m^3)
#' @return `cpw` Specific heat capacity of woody vegetation (J / kg / K)
#' @return `phw` Density of woody vegetation (no air) (kg / m^3)
#' @return `kwood` # thermal conductivity of wood (W / m / K)
#' @details if `tme` is a single value `PAI` and `pLAI` are vectors of length `m` giving
#' Plant Area Index values and proportions of green vegetation for each canopy node `m`.
#' If `tme` is a vector of times, `PAI` and `pLAI` are arrays of dimension m x number of
#' hours. I.e. seperate `PAI` and `pLAI` values are derived for each hour and node. If `tme`
#' is not in hourly time increments, it is conveted to hourly time increments, from 00:00 hrs
#' on the first day toup to and including 23:00 hrs on the final day. If `PAIt` is provided,
#' and a single value, then the length of `tme` is ignored and returned `PAI` is a vector of values for
#' each node. If `PAIt` is a vector of values of length `m` then seasonally adjusted values for
#' each node and hour of `tme` are returned. If `PAIt` is a vector of values of any length
#' other than `m` spline interpolated values for each node and hour of `tme` are returned.
#' If `PAIt` is NA (the default) PAI is estimated from habitat and location using [PAIfromhabitat()].
#' @export
habitatvars <- function(habitat, lat, long, tme, m = 20, PAIt = NA) {
  # By habitat type
  if (habitat == "Evergreen needleleaf forest" | habitat == 1) {
    m2 <- round((1/15)*m,0)
    pl<-.habgen(habitat,lat,long,tme,m,m2,7.5,70,0.75)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.4,0.7,0.1)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.25 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.3 #  Reflectivity of leaves to PAR
    lw <- 0.01
    gsmax <- 0.33
    phw <- 500
    uhgt <- 1
    wgt <- (m2 / m) * 0.25
    if (class(PAI) == "matrix") {
      sPAI<-apply(PAI,1,sum)*wgt
    } else sPAI <- sum(PAI)*wgt
    zm0 <- roughlength(uhgt, PAI = sPAI)
  }
  if (habitat == "Evergreen Broadleaf forest" | habitat == 2) {
    m2 <- round((4/20)*m,0)
    pl<-.habgen(habitat,lat,long,tme,m,m2,6.5,70,0.8)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.55, 1,0.2)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.3 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.35
    gsmax <- 0.33
    phw <- 1100
    uhgt <- 4
    wgt <- (m2 / m) * 0.25
    if (class(PAI) == "matrix") {
      sPAI<-apply(PAI,1,sum)*wgt
    } else sPAI <- sum(PAI)*wgt
    zm0 <- roughlength(uhgt, PAI = sPAI)
  }
  if (habitat == "Deciduous needleleaf forest" | habitat == 3) {
    m2 <- round((1/10)*m,0)
    pl<-.habgen(habitat,lat,long,tme,m,m2,7.5,70,0.75)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.4,0.7,0.1)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.3 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.25 #  Reflectivity of leaves to PAR
    lw <- 0.01
    gsmax <- 0.33
    phw <- 500
    uhgt <- 1
    wgt <- (m2 / m) * 0.25
    if (class(PAI) == "matrix") {
      sPAI<-apply(PAI,1,sum)*wgt
    } else sPAI <- sum(PAI)*wgt
    zm0 <- roughlength(uhgt, PAI = sPAI)
  }
  if (habitat == "Deciduous broadleaf forest" | habitat == 4) {
    m2 <- round((2/15)*m,0)
    pl<-.habgen(habitat,lat,long,tme,m,m2,6.5,70,0.775)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.5, 0.6,0.15)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.3 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.07
    gsmax <- 0.23
    phw <- 700
    uhgt <- 2
    wgt <- (m2 / m) * 0.25
    if (class(PAI) == "matrix") {
      sPAI<-apply(PAI,1,sum)*wgt
    } else sPAI <- sum(PAI)*wgt
    zm0 <- roughlength(uhgt, PAI = sPAI)
  }
  if (habitat == "Mixed forest" | habitat == 5) {
    m2 <- round((1.5/10)*m,0)
    pl<-.habgen(habitat,lat,long,tme,m,m2,7,70,0.775)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.45,0.65,0.12)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.28 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.25 #  Reflectivity of leaves to PAR
    lw <- 0.04
    gsmax <- 0.28
    phw <- 600
    uhgt <- 1.5
    wgt <- (m2 / m) * 0.25
    if (class(PAI) == "matrix") {
      sPAI<-apply(PAI,1,sum)*wgt
    } else sPAI <- sum(PAI)*wgt
    zm0 <- roughlength(uhgt, PAI = sPAI)
  }
  if (habitat == "Closed shrublands" | habitat == 6) {
    pl<-.habgen(habitat,lat,long,tme,m,1,6,80,0.85,under=F)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.2,0.5,0.05)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.3 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.03
    gsmax <- 0.35
    phw <- 500
    uhgt <- 0.05
    zm0 <- 0.004
  }
  if (habitat == "Open shrublands" | habitat == 7) {
    pl<-.habgen(habitat,lat,long,tme,m,1,6,80,0.75,under=F)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.2,0.5,0.05)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.3 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.03
    gsmax <- 0.35
    phw <- 500
    uhgt <- 0.05
    zm0 <- 0.004
  }
  if (habitat == "Woody savannas" | habitat == 8) {
    m2<-round((0.75/3)*m,0)
    pl<-.habgen(habitat,lat,long,tme,m,1,6.5,70,0.75)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.2,0.3,0.1)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.35 # reflectivity (shortwave radiation) of leaves
    refw = 0.2 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.01
    gsmax <- 0.33
    phw <- 300
    uhgt <- 0.75
    wgt <- (m2 / m) * 0.25
    if (class(PAI) == "matrix") {
      sPAI<-apply(PAI,1,sum)*wgt
    } else sPAI <- sum(PAI)*wgt
    zm0 <- roughlength(uhgt, PAI = sPAI)
  }
  if (habitat == "Savannas" | habitat == 9) {
    pl<-.habgen(habitat,lat,long,tme,m,1,1,50,0.7,under=F)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.2,0.3,0.02)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.35 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.01
    gsmax <- 0.33
    phw <- 300
    uhgt <- 0.05
    zm0 <- 0.004
  }
  if (habitat == "Short grasslands" | habitat == 10) {
    pl<-.habgen(habitat,lat,long,tme,m,1,1,50,0.85,under=F)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.2,0.3,0.02)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.35 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.01
    gsmax <- 0.33
    phw <- 300
    uhgt <- 0.05
    zm0 <- 0.004
  }
  if (habitat == "Tall grasslands" | habitat == 11) {
    pl<-.habgen(habitat,lat,long,tme,m,1,1,50,0.85,under=F)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.2,0.3,0.02)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.35 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.01
    gsmax <- 0.33
    phw <- 300
    uhgt <- 0.05
    zm0 <- 0.004
  }
  if (habitat == "Permanent wetlands" | habitat == 12) {
    pl<-.habgen(habitat,lat,long,tme,m,1,1,50,0.95,under=F)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.2,0.3,0.02)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.5 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.09
    gsmax <- 0.55
    phw <- 600
    uhgt <- 0.02
    zm0 <- 0.002
  }
  if (habitat == "Croplands" | habitat == 13) {
    pl<-.habgen(habitat,lat,long,tme,m,1,1,50,0.8,under=F)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.25,0.3,0.02)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.3 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.02
    gsmax <- 0.33
    phw <- 300
    uhgt <- 0.05
    zm0 <- 0.004
  }
  if (habitat == "Urban and built-up" | habitat == 14) {
    m2 <- round((0.1/0.8) * m, 0)
    pl<-.habgen(habitat,lat,long,tme,m,m2,6.5,70,0.75)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.5, 0.6,0.15)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.3 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.04
    gsmax <- 0.33
    phw <- 600
    uhgt <- 0.1
    wgt <- (m2 / m) * 0.25
    if (class(PAI) == "matrix") {
      sPAI<-apply(PAI,1,sum)*wgt
    } else sPAI <- sum(PAI)*wgt
    zm0 <- roughlength(uhgt, PAI = sPAI)
  }
  if (habitat == "Cropland/Natural vegetation mosaic" | habitat == 15) {
    pl<-.habgen(habitat,lat,long,tme,m,1,1,50,0.75,under=F)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.35,0.5,0.02)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.3 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.03
    gsmax <- 0.3
    phw <- 400
    uhgt <- 0.05
    zm0 <- 0.004
  }
  if (habitat == "Barren or sparsely vegetated" | habitat == 16) {
    pl<-.habgen(habitat,lat,long,tme,m,1,1,50,0.55,under=F)
    PAI<-pl$PAI
    pLAI<-pl$pLAI
    thickw <- thickgeometry(m, 0.25,0.3,0.02)
    iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
    refls = 0.3 # reflectivity (shortwave radiation) of leaves
    refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
    reflp = 0.2 #  Reflectivity of leaves to PAR
    lw <- 0.015
    gsmax <- 0.33
    phw <- 600
    uhgt <- 0.01
    zm0 <- 0.001
  }
  return(list(hgt = pai$height, PAI = PAI, x = pai$x, lw = lw, cd = 0.2, iw = iw,
              hgtg = uhgt, zm0 = zm0, pLAI = pLAI, refls = refls,
              refg = 0.15, refw = refw, reflp = reflp, vegem = 0.97, gsmax = gsmax,
              q50 = 100, thickw = thickw, cpw = 1200, phw = phw, kwood = 0.14))
}
