#' A data frame of hourly weather
#'
#' A data frame of hourly weather in 1995 at Camborne weather station, Cornwall (50.2178N, 5.32656W)
#'
#' @format a data frame with the following elements:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times}
#'  \item{temp}{temperature (degrees C)}
#'  \item{relhum}{relative humidity (percentage)}
#'  \item{pres}{atmospheric press (kPa)}
#'  \item{swrad}{Total incoming shortwave radiation (W / m^2)}
#'  \item{difrad}{Diffuse radiation (W / m^2)}
#'  \item{skyem}{Sky emissivity (0-1)}
#'  \item{windspeed}{Wind speed (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#' }
"climdata"
#' A dataset of global climate variables
#'
#' A global dataset containing containing the following climate variables averaged
#' over the period 2008 to 2017 as used by [PAIfromhabitat()].
#'
#' @format An array with 94 rows, 192 columns and the following five climate variables:
#' \describe{
#'   \item{1}{mean annual temperature (ºC)}
#'   \item{2}{coefficient of variation in temperature (K)}
#'   \item{3}{mean annual temperature (mm per year)}
#'   \item{4}{coefficient of variation in annual rainfall (mm per 0.25 days)}
#'   \item{5}{numeric month with the most rainfall (1-12)}
#' }
#' @source \url{http://www.ncep.noaa.gov/}
"globalclimate"
#' A table of habitat types
#'
#' A dataset containing habitat descriptors and habitat numbers as used by
#' [PAIfromhabitat()].
#'
#' @format a data.frame with 17 rows and two columns:
#' \describe{
#'   \item{number}{Numeric values for each habitat type}
#'   \item{descriptor}{A description of the habitat type}
#' }
#' #' @source \url{https://modis.gsfc.nasa.gov/data/}
"habitats"
