% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplemodel.R
\name{Tground}
\alias{Tground}
\title{Compute ground surface temperature}
\usage{
Tground(
  climdata,
  lat,
  long,
  hgt,
  pai,
  x = 1,
  gsmax = 0.33,
  slope = 0,
  aspect = 0,
  alb = 0.23,
  galb = 0.15,
  zu = 2,
  Bowen = 1.5,
  upper = 25,
  umin = 0.5,
  merid = 0,
  dst = 0,
  clump = 0,
  method = "S",
  lims = TRUE
)
}
\arguments{
\item{climdata}{a dataframe of hourly weather data formated and with units
as per the internal dataset \code{climdata}}

\item{lat}{the latitude of the location (decimal degrees)}

\item{long}{the longitude of the location (decimal degrees)}

\item{hgt}{canopy height (m)}

\item{pai}{the total one sided area of canopy elements per unit ground area (see details)}

\item{x}{leaf distribution angle coefficient}

\item{gsmax}{maximum stomatal conductance of leaves (mol / m^2 / s)}

\item{slope}{the slope of the underlying ground surface (decimal degrees). Ignored if \code{method} = S.}

\item{aspect}{the aspect of the underlying ground surface (decimal degrees). Ignored if \code{method} = S.}

\item{alb}{the albedo of the canopy surface (either the combined ground and
canopy albedo if \code{method = S} or just the canopy albedo)}

\item{galb}{ground surface albedo}

\item{zu}{the height above ground of wind speeds in \code{climdata} (m)}

\item{Bowen}{Optional parameter specifying the Bowen Ratio of the surface. Used to improve estimates
of temperature in application of the Penman-Monteith equation}

\item{upper}{optional upper limit to temperature offset (difference between reference
and canopy surface temperature cannot exceed this value). Ignored if \code{lims} = FALSE}

\item{umin}{optional minimum wind speed for computing conductances (avoids conductances being too low)}

\item{merid}{an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT)}

\item{dst}{an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).}

\item{clump}{clumpiness factor (0-1, see details)}

\item{method}{if \code{S} treats the vegetation surface as a flat surface (see details)}

\item{lims}{optional logical indicating whether to limit temperatures by \code{upper} and dewpoint temperature}
}
\value{
a vector of ground surface temperatures (deg C)
}
\description{
Function to compute the ground surface temperature using the Penman-Monteith equation
}
\details{
if \code{pai} is unknown it can be estimated as -ln(1-fractional
canopy cover). if \code{clump} = 0 the canopy is assumed entirely uniform
and radiation transmission is as for a turbid medium. As \code{clump}
approaches 1, the canopy is assumed to be increasingly patchy, such
that a greater proportion of reaches the ground without being obscured
by leaves. If \code{method} = S, the radiation intercepted by the canopy is assumed to
be that for a flat surface. If \code{method} is not S, the radiation absorption
by the ground and canopy surface are computed separately accounting for the
inclination of the ground surface and the distribution of leaf angles.
}
\examples{
# Compute temperature a surface
Tg<-Tground(climdata,50.2178,-5.32656,hgt=0.5,pai=2,x=1,gsmax=0.33,alb=0.23)
plot(Tg,type="l")
}
