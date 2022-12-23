#' Calculates the astronomical Julian day
#'
#' @description `jday` is used to calculate the astronomical Julian day (days since since January 1, 4713 BCE at noon UTC) from a given year, month and day.
#'
#' @param year year (AD).
#' @param month month in numeric form (1-12).
#' @param day days of the month (1-31).
#' @param hour hours (decimal, 0-23).
#' @param min minutes (decimal, 0-59).
#' @param sec seconds (decimal, 0-59).
#' @param dst an optional numeric value specifying the time zone expressed as hours different from GMT (-ve to west).
#' @param tme POSIXlt object of times in UTC. Other inputs ignored if this is provided
#'
#' @return Julian Day. I.e. the number of days since January 1, 4713 BCE at noon UTC.
#' @export
#'
#' @examples
#' tme <- as.POSIXlt(0, origin = "2020-03-21 12:43", tz = "UTC")
#' jday(tme = tme)
#' jd1 <- jday(2010, 1, 31)
#' jd2 <- jday(2010, 1, 31, 11, 0, 0)
#' jd1 - jd2
jday <- function(year, month, day, hour = 12, min = 0, sec = 0, dst = 0, tme = NA) {
  if (class(tme)[1] != "logical") {
    year<- tme$year + 1900
    month <- tme$mon + 1
    day <- tme$mday
    hour <- tme$hour
    min <- tme$min
    sec <- tme$sec
    dst <- 0
  }
  day_decimal <- day + (hour - dst + (min + sec / 60) / 60) / 24
  monthadj <- month + (month < 3) * 12
  yearadj <- year + (month < 3) * -1
  julian_day <- trunc(365.25 * (yearadj + 4716)) + trunc(30.6001 *
                                                           (monthadj + 1)) + day_decimal - 1524.5
  B <- (2 - trunc(yearadj / 100) + trunc(trunc(yearadj / 100) / 4))
  julian_day <- julian_day + (julian_day > 2299160) * B
  julian_day
}
#' Calculates the solar time
#'
#' @description `solartime` is used to calculate the solar time. I.e. the time that would be measured by a sundial.
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param long longitude of the location for which the solar time is required (decimal degrees, -ve west of Greenwich meridian).
#' @param jd Julian day expressed as an integer as returned by [jday()].
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#'
#' @return the solar time. I.e. the times that would be measured by a sundial (hours).
#' @export
#'
#' @details
#' ‘solartime’ accounts for two factors: firstly, east or west component of the analemma, namely
#' the angular offset of the Sun from its mean position on the celestial sphere as viewed from
#' Earth due the eccentricity of the Earth's orbit and the obliquity due to tilt of the Earth's
#' rotational axis. These two factors have different wavelengths, amplitudes and phases,
#' that vary over geological timescales. The equations used here are those derived by Milne. BY default,
#' local are used, with the meridian set to round(long / 15, 0) * 15.
#'
#'
#' @examples
#' jd <- jday (2010, 6, 21) # Julian day
#' solartime(12, -5, jd) # solartime at noon on 21 June 2010, 5ºW
solartime <- function(localtime, long, jd, merid = round(long / 15, 0) * 15, dst = 0) {
  m <- 6.24004077 + 0.01720197 * (jd - 2451545)
  eot <- -7.659 * sin(m) + 9.863 * sin (2 * m + 3.5932)
  st <- localtime + (4 * (long - merid) + eot) / 60 - dst
  st
}
#' Calculates the solar azimuth
#'
#' @description `solazi` is used to calculate the solar azimuth at any given location from the local time.
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param lat latitude of the location for which the solar azimuth is required (decimal degrees, -ve south of the equator).
#' @param long longitude of the location for which the solar azimuth is required (decimal degrees, -ve west of Greenwich meridian).
#' @param jd Julian day expressed as an integer as returned by [jday()].
#' @param merid optional  value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst optional value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @return a numeric value or vector of values representing the solar azimuth (decimal degrees).
#' @export
#'
#' @examples
#' # solar azimuth at noon on 21 June 2010, Porthleven, Cornwall, UK
#' jd <- jday (2010, 6, 21) # Julian day
#' solazi(12, 50.08, -5.31, jd)
solazi <- function(localtime, lat, long, jd, merid = round(long / 15, 0) * 15, dst = 0) {
  stime <- solartime(localtime, long, jd, merid, dst)
  tt <- 0.261799 * (stime - 12)
  declin <- (pi * 23.5 / 180) * cos(2 * pi * ((julian - 171) / 365.25))
  Sinh <- sin(declin) * sin(lat * pi / 180) + cos(declin) *
    cos(lat * pi / 180) * cos(tt)
  hh <- (atan(Sinh / sqrt(1 - Sinh * Sinh)))
  sinazi <- cos(declin) * sin(tt) / cos(hh)
  cosazi<-(sin(lat*pi/180)*cos(declin)*cos(tt)-cos(pi*lat/180)*sin(declin))/
           sqrt((cos(declin)*sin(tt))^2+(sin(pi*lat/180)*cos(declin)*cos(tt)-cos(pi*lat/180) * sin(declin)) ^ 2)
  sqt <- 1 - sinazi * sinazi
  sqt[sqt < 0] <- 0
  solz <- 180 + (180 * atan(sinazi / sqrt(sqt))) / pi
  solz[cosazi < 0 & sinazi < 0] <- 180 - solz[cosazi < 0 & sinazi < 0]
  solz[cosazi < 0 & sinazi >= 0] <- 540 - solz[cosazi < 0 & sinazi >= 0]
  solz
}
#' Calculates the solar altitude
#'
#' @description `solalt` is used to calculate the solar altitude at any given location from the local time.
#'
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param lat latitude of the location for which the solar altitude is required (decimal degrees, -ve south of the equator).
#' @param long longitude of the location for which the solar altitude is required (decimal degrees, -ve west of Greenwich meridian).
#' @param jd Julian day expressed as an integer as returned by [jday()].
#' @param merid optional value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst optional value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#'
#' @return solar altitude (decimal º).
#' @export
#'
#' @examples
#' # solar altitude at noon on 21 June 2010, Porthleven, Cornwall
#' jd <- jday (2010, 6, 21) # Julian day
#' solalt(12, 50.08, -5.31, jd)
solalt <- function(localtime, lat, long, jd, merid = round(long / 15, 0) * 15, dst = 0) {
  stime <- solartime(localtime, long, jd, merid, dst)
  tt <- 0.261799 * (stime - 12)
  declin <- (pi * 23.5 / 180) * cos(2 * pi * ((julian - 171) / 365.25))
  sinh <- sin(declin) * sin(lat * pi / 180) + cos(declin) *
    cos(lat * pi / 180) * cos(tt)
  sa <- (180 * atan(sinh / sqrt(1 - sinh * sinh))) / pi
  sa
}
#' Calculates the solar coefficient
#'
#' @description `solarcoef` is used to calculate the proportion of direct beam radiation
#' incident on an inclined surface at a specified time and location.
#'
#' @param slope slopes angle (decimal º)
#' @param aspect aspect angle (decimal º)
#' @param localtime local time (decimal hour, 24 hour clock).
#' @param lat latitude of the location for which the coefficientis required (decimal degrees, -ve south of the equator).
#' @param long longitude of the location for which the solar index is required (decimal degrees, -ve west of Greenwich meridian).
#' @param jd integer representing the Julian day as returned by [jday()].
#' @param merid optional value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param horizon optinal logical indicating whether to set coefficient to zero if solar altitude < 0.
#' @return the proportion of direct beam radiation incident on an inclined surface.
#' @export
#'
#' @seealso the microclimafunction [microclima::solarindex()] can be used to derive the solar
#' coeffficient for cells of a digital elevation dataset accounting for topographic shading.
#'
#' @examples
#' library(raster)
#' jd <- jday (2010, 6, 21) # Julian day
#' lt <- c(0:2400) / 100
#' si <- solarcoef(0, 0, lt, 50, -5, jd)
#' plot(si ~ lt, type = "l")
solarcoef <- function(slope, aspect, localtime, lat, long, jd,
                      merid = round(long / 15, 0) * 15, dst = 0,
                      horizon = TRUE) {
  saltitude <- solalt(localtime, lat, long, jd, merid, dst)
  alt <- saltitude * (pi / 180)
  zen <- pi / 2 - alt
  sazimuth <- solazi(localtime, lat, long, jd, merid, dst)
  azi <- sazimuth * (pi / 180)
  sl <- slope * (pi / 180)
  asp <- aspect * (pi / 180)
  index <- cos(zen) * cos(sl) + sin(zen) * sin(sl) * cos(azi - asp)
  index[index < 0] <- 0
  index[alt < 0] <- 0
  index
}
#' Calculates clear sky radiation
#'
#' @param tme POSIXlt object of times in UTC.
#' @param lat latitude of the location for which the solar azimuth is required (decimal degrees, -ve south of the equator).
#' @param long longitude of the location for which the solar azimuth is required (decimal degrees, -ve west of Greenwich meridian).
#' @param h optional specific humidities (\ifelse{html}{\out{kg kg<sup>{-1}</sup> }}{\eqn{kg kg^{-1}}}).
#' @param tc optional temperatures (ºC).
#' @param pk optional pressure (kPa).
#' @param G optional value describing he moisture profile in the atmosphere (per Smith 1966).
#' @param Ie an optional single value for extra-terrestrail radiation to permit adjustment for
#' sun-earth distances
#' @param merid optional  value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst optional value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @return expected clear-sky radioation (W/m^2.
#' @export
clearskyrad <- function(tme, lat, long, h = 0.00697, tc = 15, pk = 101.3, G = 2.78, Ie = 1352.778,
                        merid = 0, dst = 0) {

  p <- pk*1000
  jd<-jday(tme=tme)
  lt <- tme$hour + tme$min / 60 + tme$sec / 3600
  sa <- solalt(lt, lat, long, jd, merid, dst)
  sa[sa < 0] <- NA
  z <- (90 - sa) * (pi / 180)
  m <- 35 * cos(z) * ((1224 * cos(z)^2 + 1)^(-0.5))
  TrTpg <- 1.021 - 0.084 * (m * 0.000949 * 0.01 * p + 0.051)^0.5
  pk <- p / 1000
  e0 <- 0.6108 * exp(17.27 * tc/(tc + 237.3))
  ws <- 0.622 * e0/pk
  rh <- (h/ws) * 100
  rh <- ifelse(rh > 100, 100, rh)
  xx <- log(rh / 100) + ((17.27 * tc) / (237.3 + tc))
  Td <- (237.3 * xx) / (17.27 - xx)
  u <- exp(0.1133 - log(G + 1) + 0.0393 * Td)
  Tw <- 1 - 0.077 * (u * m) ^ 0.3
  Ta <- 0.935 * m
  od <- TrTpg * Tw * Ta
  Ic <- Ie * (cos(z)) * od
  Ic[Ic > Ie] <- NA
  Ic[Ic < 0] <- NA
  Ic
}
#' Calculates direct beam radiation transmission through vegetated canopies
#'
#' @param l leaf area index (total one-sided leaf area per unit ground area)
#' @param x leaf distribution angle coefficient
#' @param sa solar altitude in decimal degrees as returned by [solalt()].
#' @param ref average reflectivity of leaves, usually in shortwave spectrum.
#' @param clump clumpiness factor (0-1, see details)
#' @return the proportion of direct beam radiation transmitted through the canopy
#' @export
#' @details if `clump` = 0 the canopy is assumed entirely uniform and radiation
#' transmission is as for a turbid medium. As `clump` approaches 1, the canopy
#' is assumed to be increasingly patchy, such that a greater proportion of reaches
#' the ground without being obscured by leaves.
#' @seealso [cansw()] to calculate total shortwave radiation underneath canopies
cantransdir <- function(l, x, sa, ref = 0.5, clump = 0) {
  f <- 1 / (1 - clump)
  zen <- 90 - sa
  k <- sqrt((x^2 + (tan(zen * (pi/180))^2)))/(x + 1.774 * (x + 1.182)^(-0.733))
  s <- sqrt(1 - ref)
  ks <- k * s
  tr <- exp(-ks * l * f)
  tr <- (1 - clump) * tr + clump
  tr
}
#' Calculates diffuse radiation transmission through vegetated canopies
#'
#' @param l leaf area index (total one-sided leaf area per unit ground area)
#' @param ref average reflectivity of leaves, usually in shortwave spectrum.
#' @param clump clumpiness factor (0-1, see details)
#' @return the proportion of diffuse beam radiation transmitted through the canopy
#' @export
#' @details if `clump` = 0 the canopy is assumed entirely uniform and radiation
#' transmission is as for a turbid medium. As `clump` approaches 1, the canopy
#' is assumed to be increasingly patchy, such that a greater proportion of reaches
#' the ground without being obscured by leaves.
#' @seealso [cansw()] to calculate total shortwave radiation underneath canopies
cantransdif <- function(l, ref = 0.25, clump = 0) {
  f <- 1 / (1 - clump)
  s <- sqrt(1 - ref)
  tr <- exp(-s * l * f)
  tr <- (1 - clump) * tr + clump
  tr
}
#' Calculates the diffuse fraction from incoming shortwave radiation
#'
#' @description `difprop` calculates proportion of incoming shortwave radiation that is diffuse radiation using the method of Skartveit et al. (1998) Solar Energy, 63: 173-183.
#'
#' @param rad a vector of incoming shortwave radiation values (either \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} or \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}})
#' @param jd the Julian day as returned by [jday()]
#' @param localtime a single numeric value representing local time (decimal hour, 24 hour clock)
#' @param lat a single numeric value representing the latitude of the location for which partitioned radiation is required (decimal degrees, -ve south of equator).
#' @param long a single numeric value representing the longitude of the location for which partitioned radiation is required (decimal degrees, -ve west of Greenwich meridian).
#' @param hourly specifies whether values of `rad` are hourly (see details).
#' @param watts a logical value indicating  whether the units of `rad` are \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}} (TRUE) or \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} (FALSE).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param corr an optional numeric value representing a correction to account for over- or under-estimated diffuse proportions. Values > 1 will apportion a greater ammount of total radiation as diffuse than originally calculated by the formula.
#'
#' @return a vector of diffuse fractions (either \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} or \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}}).
#' @export
#'
#' @details
#' The method assumes the environment is snow free. Both overall cloud cover and heterogeneity in
#' cloud cover affect the diffuse fraction. Breaks in an extensive cloud deck may primarily
#' enhance the beam irradiance, whereas scattered clouds may enhance the diffuse irradiance and
#' leave the beam irradiance unaffected.  In consequence, if hourly data are available, an index
#' is applied to detect the presence of such variable/inhomogeneous clouds, based on variability
#' in radiation for each hour in question and values in the preceding and deciding hour.  If
#' hourly data are unavailable, an average variability is determined from radiation intensity.
#'
#' @examples
#' rad <- c(5:42) / 0.036 # typical values of radiation in W/m^2
#' jd <- jday(2017, 6, 21) # julian day
#' dfr <- difprop(rad, jd, 12, 50, -5)
#' plot(dfr ~ rad, type = "l", lwd = 2, xlab = "Incoming shortwave radiation",
#'      ylab = "Diffuse fraction")
difprop <- function(rad, jd, localtime, lat, long, hourly = FALSE,
                    watts = TRUE, merid = round(long / 15, 0) * 15, dst = 0,
                    corr = 1) {
  if (watts) rad <- rad * 0.0036
  sa <- solalt(localtime, lat, long, jd, merid, dst)
  alt <- sa * (pi / 180)
  k1 <- 0.83 - 0.56 * exp(- 0.06 * sa)
  si <- cos(pi / 2 - alt)
  si[si < 0] <- 0
  k <- rad / (4.87 * si)
  k[!is.finite(k)] <- 0
  k <- ifelse(k > k1, k1, k)
  k[k < 0] <- 0
  rho <- k / k1
  if (hourly) {
    rho <- c(rho[1], rho, rho[length(rho)])
    sigma3  <- 0
    for (i in 1:length(rad)) {
      sigma3[i] <- (((rho[i + 1] - rho[i]) ^ 2 + (rho[i + 1] - rho[i + 2]) ^ 2)
                    / 2) ^ 0.5
    }
  } else {
    sigma3a <- 0.021 + 0.397 * rho - 0.231 * rho ^ 2 - 0.13 *
      exp(-1 * (((rho - 0.931) / 0.134) ^ 2) ^ 0.834)
    sigma3b <- 0.12 + 0.65 * (rho - 1.04)
    sigma3 <- ifelse(rho <= 1.04, sigma3a, sigma3b)
  }
  k2 <- 0.95 * k1
  d1 <- ifelse(sa > 1.4, 0.07 + 0.046 * (90 - sa) / (sa + 3), 1)
  K <- 0.5 * (1 + sin(pi * (k - 0.22) / (k1 - 0.22) - pi / 2))
  d2 <- 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 * K ^ 2))
  d3 <- (d2 * k2) * (1 - k) / (k * (1 - k2))
  alpha <- (1 / sin(alt)) ^ 0.6
  kbmax <- 0.81 ^ alpha
  kmax <- (kbmax + d2 * k2 / (1 - k2)) / (1 + d2 * k2 / (1 - k2))
  dmax <- (d2 * k2) * (1 - kmax) / (kmax * (1 - k2))
  d4 <- 1 - kmax * (1 - dmax) / k
  d <- ifelse(k <= kmax, d3, d4)
  d <- ifelse(k <= k2, d2, d)
  d <- ifelse(k <= 0.22, 1, d)
  kX <- 0.56 - 0.32 * exp(-0.06 * sa)
  kL <- (k - 0.14) / (kX - 0.14)
  kR <- (k - kX) / 0.71
  delta <- ifelse(k >= 0.14 & k < kX, -3 * kL ^ 2 *(1 - kL) * sigma3 ^ 1.3, 0)
  delta <- ifelse(k >= kX & k < (kX + 0.71), 3 * kR * (1 - kR) ^ 2 * sigma3 ^
                    0.6, delta)
  d[sigma3 > 0.01] <- d[sigma3 > 0.01] + delta[sigma3 > 0.01]
  d[rad == 0] <- 0.5
  d[sa < 0] <- 1
  # apply correction
  dif_val <- rad * d
  dif_val_adj <- dif_val * corr
  d <- dif_val_adj /rad
  d[d > 1] <- 1
  d[d < 0] <- 1
  d[!is.finite(d)] <- 0.5
  d
}
#' Calculates radiation received under vegetated canopies
#'
#' @description `cansw` calculates the flux density of incoming shortwave radiation underneath vegetation canopies
#'
#' @param globrad a vector of incoming shortwave radiation values (either \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} or \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}})
#' @param dp optional numeric value or vector of values indicating the proportion of `globarad` that is diffuse.
#' If no value is provided, calculated using [difprop()]
#' @param jd the Julian day as returned by [jday()]
#' @param localtime a single numeric value representing local time (decimal hour, 24 hour clock)
#' @param lat latitude of location (decimal degrees, -ve south of equator).
#' @param long longitude of location (decimal degrees, -ve west of Greenwich meridian).
#' @param l leaf area index (total one-sided area of leaf per unit ground area)
#' @param x leaf distribution angle coefficient
#' @param ref average reflectivity of leaves in shortwave spectrum.
#' @param hourly Used if dp = NA. Specifies whether values of `rad` are hourly (see details).
#' @param watts Used if dp = NA. A logical value indicating  whether the units of `rad` are \ifelse{html}{\out{W m<sup>-2</sup>}}{\eqn{W m^{-2}}} (TRUE) or \ifelse{html}{\out{MJ m<sup>-2</sup> hr<sup>-1</sup>}}{\eqn{MJ m^{-2} hr^{-1}}} (FALSE).
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT). Default is `round(long / 15, 0) * 15`
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param corr Used if dp = NA. An optional numeric value representing a correction to account for over- or under-estimated diffuse proportions. Values > 1 will apportion a greater ammount of total radiation as diffuse than originally calculated by the formula.
#' @param tme object of POSIXlt indicating times in UTC. Can be used in place of specifying
#' `jd`, `localtime`, `dst` and `merid`
#' @param clump clumpiness factor for canopy (0-1, see details)
#' @return shortwave radiation received underneath vegetated canopies. Units as for globrad.
#' @export
#' @details Calculated the flux density of radiation received below leaf area l, not absorbed.
#' I.e. ground albedo not accounted for. if `clump` = 0 the canopy is assumed entirely
#' uniform and radiation transmission is as for a turbid medium. As `clump` approaches 1,
#' the canopy is assumed to be increasingly patchy, such that a greater proportion of reaches
#' the ground without being obscured by leaves.
#'
#' @examples
#' l <- c(0:1000)/100
#' sw <- cansw(500, dp = 0.5, jday(2020, 6, 21), 12, 50, -5, l, 1)
#' plot(sw ~ l, type = "l")
#' @seealso [cantransdir()], [cantransdif()] and [canlw()]
#'
cansw <- function(globrad, dp = NA, jd, localtime, lat, long, l, x, ref = 0.2,
                  hourly = FALSE, watts = TRUE, merid = round(long / 15, 0) * 15,
                  dst = 0, corr = 1, tme = NA, clump = 0) {
  if (class(tme)[1] != "logical") {
    jd <- jday(tme = tme)
    localtime <- tme$hour + tme$min / 60 + tme$sec / 3600
    dst <- 0
    merid <- 0
  }
  if (class(dp) == "logical") {
    dp <- difprop(globrad, jd, localtime, lat, long, hourly, watts, merid, corr)
  }
  sa <- solalt(localtime, lat, long, jd, merid, dst)
  drtr <- cantransdir(l, x, sa, ref, clump)
  ditr <- cantransdif(l, ref, clump)
  rad <- dp * globrad * ditr + (1 - dp) * globrad * drtr
  rad
}
#' Calculates longwave radiation underneath vegetated canopies
#'
#' @description `canlw` calculates the flux density (W / m^2) of incoming and
#' outgoing longwave radiation underneath vegetated canopies
#'
#' @param tc temperature (deg C)
#' @param l leaf area index
#' @param ref average leaf reflectivity in longwave spectrum (1 - thermal emissivity)
#' @param skyem sky emissivity
#' @param clump clumpiness factor for canopy (0-1, see details)
#' @return A list wiht the following elements:
#' @return `lwout` Flux density of total outgoing longwave radiation emitted under leaf area `l` (W / m2)
#' @return `lwin` Flux density of incoming longwave radiation received under leaf area `l` (W / m2)
#' @return `lwabs` Flux density of absorbed longwave radiation under leaf area `l` (W / m2)
#' @return `lwnet` Flux density of net longwave radiation emitted under leaf area `l` (W / m2)
#' @export
#' @examples
#' l <- c(0:1000) / 500
#' lw1 <- canlw(11, l, skyem = 0.9)
#' lw2 <- canlw(11, l, skyem = 0.7)
#' lw3 <- canlw(11, l, skyem = 0.5)
#' plot(lw1$lwnet ~ l, type = "l", lwd = 2, ylim = c(0,200), ylab = "Net longwave")
#' par(new = TRUE)
#' plot(lw2$lwnet ~ l, type = "l", col = "blue", lwd = 2, ylim = c(0,200), ylab = "")
#' par(new = TRUE)
#' plot(lw3$lwnet ~ l, type = "l", col = "red", lwd = 2, ylim = c(0,200), ylab = "")
canlw <- function(tc, l, ref = 0.03, skyem = 0.9, clump = 0) {
  f <- 1 / (1 - clump)
  lwout <- (1 - ref) * 5.67 * 10^-8 * (tc + 273.15) ^ 4
  s <- sqrt(1 - ref)
  tr <- exp(-s * l * f)
  tr <- (1 - clump) * tr + clump
  lwcan <- (1 - tr) * lwout
  lwsky <- skyem * tr^2 * lwout
  lwin <- lwcan + lwsky
  return(list(lwout = lwout, lwin = lwin,
              lwabs = (1 - ref) * lwin,
              lwnet = lwout - (1 - ref) * lwin))
}
#' Calculates a multiplication factor for radiation absorbtion by sunlit leaves
#'
#' @description `radmult` calaculates a leaf angle-dependent multiplication factor for sunlit
#' leaves.
#' @param x leaf angle distribution coefficient (ratio nof horizontal to vertical projection of foliage)
#' @param sa solar altitude (decimal degrees)
#' @return multiplication factor representing ratio of radiation intercepted by leaves relative
#' to a horizontal surface
#' @details An inclined leaf facing in the direction of the sun will tend to intercept more direct
#' beam radiation than a horizontal surface. This function uses an approximation of the
#' x-dependent leaf-angle density function as detailed in Campbell (1990) Agricultural and Forest Meteorology, 49: 173-176
#' to calculate multiplication factor representing ratio of radiation intercepted by leaves relative
#' to a horizontal surface.
#' @export
#' @seealso [psunlit()]
radmult <- function(x, sa) {
  sa<-ifelse(sa<0,0,sa)
  coef1 <- 1.20619*x^0.40711-4.88761
  coef2 <- -0.41215*x^ 0.31701+1.32412
  lrat<-coef1+coef2*log(sa)
  mult <- 1/exp(lrat)
  mult
}
#' Calculates the fraction of sunlit leaves
#'
#' @description `psunlit` calculates the fraction of sunlit leaves for application of `radmult`
#' @param l leaf area index
#' @param x leaf angle distribution coefficient (ratio nof horizontal to vertical projection of foliage)
#' @param sa solar altitude (decimal degrees)
#' @param clump clumpiness factor for canopy (0-1, see details)
#' @export
#' @seealso [radmult()]
#' @details if `clump` = 0 the canopy is assumed entirely uniform and the sunlit proportion
#' is as for a turbid medium. As `clump` approaches 1, the canopy is assumed to be
#' increasingly patchy, such that a greater proportion of the canopy is sunlit.
psunlit <- function(l, x, sa, clump = 0) {
  f <- 1 / (1 - clump)
  sa<-ifelse(sa<0,0,sa)
  ze<-90-sa
  K <- sqrt((x^2 + (tan(ze * (pi/180))^2)))/(x + 1.774 * (x + 1.182)^(-0.733))
  Ls<-(1-exp(-K*l*f))/K
  Lp <- Ls/l
  Lp <- (1 - clump) * Lp + clump
  Lp[is.na(Lp)]<-1
  Lp
}
#' Calculates the effective fraction of sunlit leaves for diffuse radiation
.psunlitd <- function(l, x, clump = 0) {
  f <- 1 / (1 - clump)
  Ls<-(1-exp(-l*f))
  Lp <- Ls/l
  Lp <- (1 - clump) * Lp + clump
  Lp[is.na(Lp)]<-1
  Lp
}
