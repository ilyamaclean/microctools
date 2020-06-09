#' Convert between different measures of humidity
#'
#' @description converts between different humidity measures, namely
#' vapour pressure or relative, absolute or specific humidity.
#'
#' @param h humidity value(s). Units as follows: specific humidity (\ifelse{html}{\out{kg kg<sup>-1</sup>}}{\eqn{kg kg^{-1}}}),
#' absolute humidity (\ifelse{html}{\out{kg m<sup>-3</sup> }}{\eqn{kg m^{-3}}}),
#' relative humidity (\%), vapour pressure (kPa).
#' @param intype a character string description of the humidity type of `h`. One of "relative", "absolute" or "specific".
#' @param tc A numeric value specifying the temperature (ºC).
#' @param pk An optional numeric value specifying the atmospheric pressure (kPa).
#'
#' @details This function converts between vapour pressure and specific,
#' relative and absolute humidity, based on sea-level pressure and
#' temperature. It returns a list of relative, absolute and specific
#' humidity and vapour pressure. If one or more of the
#' relative humidity values exceeds 100\% a warning is given.
#'
#' @return a list of numeric humidity values with the following components:
#' @return `relative` relative humidity (\%).
#' @return `absolute`  absolute humidity (\ifelse{html}{\out{kg m<sup>-3</sup> }}{\eqn{kg m^{-3}}}).
#' @return `specific` specific humidity (\ifelse{html}{\out{kg kg<sup>-1</sup> }}{\eqn{kg kg^{-1}}}).
#' @return `vapour_pressure` vapour pressure (kPa).
#' @export
#'
#' @examples
#' converthumidity(90, 'relative', 20)
#' converthumidity(0.01555486, 'absolute', 20)
#' converthumidity(0.01292172, 'specific', 20)
converthumidity <- function (h, intype = "relative", tc = 11, pk = 101.3) {
  tk <- tc + 273.15
  if (intype != "specific" & intype != "relative" & intype !=
      "absolute") {
    warning("No valid input humidity specified. Humidity assumed to be\n relative")
    intype <- "relative"
  }
  e0 <- 0.6108 * exp(17.27 * tc/(tc + 237.3))
  ws <- (18.02 / 28.97) * (e0 / pk)
  if (intype == "specific") {
    hr <- (h/ws) * 100
  }
  if (intype == "absolute") {
    ea <- (tk * h) / 2.16679
    hr <- (ea/e0) * 100
  }
  if (intype == "relative")
    hr <- h
  if (max(hr, na.rm = T) > 100)
    warning(paste("Some relative humidity values > 100%",
                  max(hr, na.rm = T)))
  if (intype == "vapour pressure") {
    hr <- (ea / e0) * 100
  }
  ea <- e0 * (hr / 100)
  hs <- (hr / 100) * ws
  ha <- 2.16679 * (ea/tk)
  return(list(relative = hr, absolute = ha, specific = hs,
              vapour_pressure = ea))
}
#' Calculate relative humidity of soil
#'
#' @description calculates relative humidity of soil from volumetric water content
#' @param theta volumetric soil moisture content (m^3 / m^3)
#' @param b soil hydraulic shape parameter
#' @param Psie soil water potential paramater (J / kg)
#' @param Smax Saturated volumetric soil moiosture content (m3 /m3)
#' @param tc soil temperature (deg C)
#' @return soil relative humdidity expressed as fraction
#' @description Default values are for Clay-loam
#' @export
soilrh <- function(theta, b = 5.2, Psie = -2.6, Smax = 0.419, tc = 11) {
  matric <- Psie*(theta/Smax)^-b
  hr<-exp((0.018*matric)/(8.31*(tc+273.15)))
  hr
}
#' Calculate canopy layer vapour conductivity
#'
#' @description adjusts stomatal conductivity based on available incoming shortwave
#' radiation
#' @param Rsw incoming shortwave radiation (W / m^2)
#' @param gsmax maximum stomatal conductivity (mol / m^2 / sec)
#' @param q50 amount of photosynthetically active radiation when stomatal conductance is at 50 percent of its maximum
#' @return stomatal conductance (mol / m^2 / sec)
#' @export
layercond <- function(Rsw, gsmax, q50 = 100) {
  rpar <- Rsw * 4.6
  gs <- (gsmax * rpar) / (rpar + q50)
  gs
}
#' Calculate saturated vapour pressure
#'
#' @description Calculates saturated vapour pressure
#' @param tc temperature (degrees C)
#' @param ice optional logical indicating whether to calculate saturated vapour pressure
#' over ice (TRUE = Yes)
#' @return saturated vapour pressure (kPa)
#' @export
#' @examples
#' # Plot saturated vapour pressure curve
#' tc <- c(-20:30)
#' es <- satvap(tc)
#' plot(es~tc, type = "l")
#' # Difference between saturated vapour pressure of water and ice
#' tc <- c(-50:0)
#' es1 <- satvap(tc)
#' es2 <- satvap(tc, ice = TRUE)
#' plot(es1 ~ tc, type = "l", lwd = 2, ylim = c(0,0.7), col = "red")
#' par(new = T)
#' plot(es2 ~ tc, type = "l", lwd = 2, ylim = c(0,0.7), col = "blue")
satvap <- function(tc, ice = FALSE) {
  if (ice) {
    e0 <- 610.78/1000
    L <- 2.834*10^6
    T0 <- 273.16
  } else {
    e0 <- 611.2/1000
    L <- (2.501*10^6) - (2340 * tc)
    T0 <- 273.15
  }
  Rv <- 461.5
  estl <- e0 * exp((L / Rv) * (1/T0 - 1/(tc + 273.15)))
  estl
}
#' Calculates dew point temperature
#'
#' @description Calculates  dew or frost point temperature
#' @param ea vapour pressure (kPa)
#' @param air temperature (degrees C)
#' @param ice optional logical indicating whether to return frost point temperature if output
#' is less than zero degrees C (TRUE = yes)
#' @return dew or frost point temperature (degrees C)
#' @examples
#' # Comparison of forst and dew point
#' ea <- c(10:100) / 50
#' tdew <- dewpoint(ea, ice = F)
#' tfrost <- dewpoint(ea)
#' plot(tfrost~ea, type = "l", col = "blue", ylim = c(-20,20))
#' par(new = T)
#' plot(tdew~ea, type = "l", col = "red", ylim = c(-20,20))
dewpoint <- function(ea, tc = 11, ice = TRUE) {
  e0 <- 611.2/1000
  L <- (2.501*10^6) - (2340 * tc)
  T0 <- 273.15
  Rv <- 461.5
  it <- 1/T0 - (Rv/L) * log(ea/e0)
  Tdew <- 1/it - 273.15
  if (ice) {
    e0 <- 610.78/1000
    L <- 2.834*10^6
    T0 <- 273.16
    it <- 1/T0 - (Rv/L) * log(ea/e0)
    Tfrost <- 1/it - 273.15
    Tdew[Tdew < 0] <- Tfrost[Tdew < 0]
  }
  Tdew
}
