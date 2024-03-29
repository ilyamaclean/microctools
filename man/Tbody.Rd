% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplemodel.R
\name{Tbody}
\alias{Tbody}
\title{Compute ectotherm body temperature}
\usage{
Tbody(
  climdata,
  z,
  dlen,
  dhgt,
  lat,
  long,
  hgt,
  pai,
  paia = NA,
  x = 1,
  gsmax = 0.33,
  lw = 0.05,
  slope = 0,
  aspect = 0,
  alb = 0.23,
  galb = 0.15,
  dalb = 0.23,
  skinwet = 1,
  zu = 2,
  Bowen = 1.5,
  upper = 25,
  umin = 0.5,
  merid = 0,
  dst = 0,
  clump = 0,
  method = "C",
  lims = TRUE
)
}
\arguments{
\item{climdata}{a dataframe of hourly weather data formated and with units
as per the internal dataset \code{climdata}}

\item{z}{the height (below canopy) for which temperature estimates are required (m)}

\item{dlen}{length of animal (m)}

\item{dhgt}{height of animal (m)}

\item{lat}{the latitude of the location (decimal degrees)}

\item{long}{the longitude of the location (decimal degrees)}

\item{hgt}{canopy height (m)}

\item{pai}{the total one sided area of canopy elements per unit ground area (see details)}

\item{paia}{optionally, the total one sided area of canopy elements per unit ground area above z (see details)}

\item{x}{leaf distribution angle coefficient}

\item{gsmax}{maximum stomatal conductance of leaves (mol / m^2 / s)}

\item{lw}{average leaf width (m)}

\item{slope}{the slope of the underlying ground surface (decimal degrees). Ignored if \code{method} = S.}

\item{aspect}{the aspect of the underlying ground surface (decimal degrees). Ignored if \code{method} = S.}

\item{alb}{the albedo of the canopy surface (either the combined ground and
canopy albedo if \code{method = S} or just the canopy albedo)}

\item{galb}{ground surface albedo}

\item{dalb}{albedo of animal}

\item{skinwet}{proportion of animal surface acting like a wet surface}

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
a vector of animal body temperatures
}
\description{
Function to compute ectotherm body temperature
}
\details{
if \code{pai} is unknown it can be estimated as -ln(1-fractional
canopy cover). If \code{paia} is unspecified vertically uniform foliage distribution is assumed.
if \code{clump} = 0 the canopy is assumed entirely uniform
and radiation transmission is as for a turbid medium. As \code{clump}
approaches 1, the canopy is assumed to be increasingly patchy, such
that a greater proportion of reaches the ground without being obscured
by leaves. If \code{method} = S, the radiation intercepted by the canopy is assumed to
be that for a flat surface. If \code{method} is not S, the radiation absorption
by the ground and canopy surface are computed separately accounting for the
inclination of the ground surface and the distribution of leaf angles.
}
\examples{
# Compute body temperatures of a toad near the surface below canopy
Ttoad<-Tbody(climdata,0.05,0.08,0.06,50.2178,-5.32656,0.5,2)
plot(Ttoad,type="l")
}
