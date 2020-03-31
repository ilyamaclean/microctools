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
#' A table of soil parameters
#'
#' A table of soil paramaters for different soil types
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Soil.type}{description of soil type}
#'   \item{Smax}{Volumetric water content at saturation (m^3 / m^3)}
#'   \item{Smin}{Residual water content (m^3 / m^3)}
#'   \item{alpha}{Shape parameter of the van Genuchten model (cm^-1)}
#'   \item{n}{Pore size distribution parameter (dimensionless, > 1)}
#'   \item{Ksat}{Saturated hydraulic conductivity (cm / day)}
#'   \item{Vq}{Volumetric quartz content of soil}
#'   \item{Vm}{Volumetric mineral content of soil}
#'   \item{Vo}{Volumetric organic content of soil}
#'   \item{Mc}{Mass fraction of clay}
#'   \item{rho}{Soil bulk density (Mg / m^3)}
#'   \item{b}{Shape parameter for Campbell model (dimensionless, > 1)}
#'   \item{psi_e}{Matric potential (J / m^3)}
#' }
#' @source: \url{https://onlinelibrary.wiley.com/doi/full/10.1002/ird.1751}
"soilparams"
#'
#' A list of climate variables
#'
#' A list of climate variables needed to run the below canopy microclimate model for one time step
#'
#' @format a list with the following elements:
#' \describe{
#'  \item{tair}{air temperature at reference height (deg C)}
#'  \item{relhum}{relative humidity at reference height (percentage)}
#'  \item{pk}{Atmospheric pressure (kPa)}
#'  \item{u}{Wind speed at reference height (m/s)}
#'  \item{tsoil}{temperature of deepest soil layer (deg C)}
#'  \item{skyem}{sky emissivity}
#'  \item{Rsw}{Incoming shortwave radiation (W/m^2)}
#'  \item{dp}{Diffuse fraction of incoming shortwave radiation}
#'  \item{psi_h}{diabatic correction factor for above canopy heat transfer}
#'  \item{psi_m}{diabatic correction factor for above canopy momentum transfer}
#'  \item{phi_m}{diabatic correction factor for below canopy momentum transfer}
#' }
"climvars"
