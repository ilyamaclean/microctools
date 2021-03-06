% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heat.R
\name{zeroplanedis}
\alias{zeroplanedis}
\title{Calculate zero plane displacement}
\usage{
zeroplanedis(hgt, PAI, cdl = 7.5)
}
\arguments{
\item{hgt}{canopy height (m)}

\item{PAI}{Plant area index}

\item{cdl}{Control parameter scaling d/h (from Raupach 1994)}
}
\value{
zero plane displacement (m)
}
\description{
Calculate zero plane displacement
}
\references{
Raupach, M.R. (1994) Simplified expressions for vegetation roughness length and zero-plane
displacement as functions of canopy height and area index. Boundary-Layer Meteorology 71: 211-216.
}
