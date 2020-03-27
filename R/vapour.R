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
#' humidityconvert(90, 'relative', 20)
#' humidityconvert(0.01555486, 'absolute', 20)
#' humidityconvert(0.01292172, 'specific', 20)
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
