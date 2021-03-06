% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{layeraverage}
\alias{layeraverage}
\title{Merges canopy air layers and calculates average properties}
\usage{
layeraverage(lmm, tc, hgt, gha, gt, zla, z, Vo, ea, X, vden, pk, PAI, TT)
}
\arguments{
\item{lmm}{list of layers to merge as returned by \code{\link[=layermerge]{layermerge()}}}

\item{tc}{air temperature of unmerged canopy layers (dec C)}

\item{hgt}{height of canopy (m)}

\item{gha}{molar heat conductance between leaf and air (mol / m^2 / sec)}

\item{gt}{molar heat conductance between air layers due to turbulent convection (mol / m^2 / sec)}

\item{zla}{mean leaf-air distance (m)}

\item{z}{height above ground of canopy nodes (m)}

\item{Vo}{mole fraction of water vapour in air in previous time-step (mol / mol)}

\item{ea}{vapur pressure in current timestep}

\item{X}{temperature to add during timestep (deg C)}

\item{vden}{Volumetric density of vegetation (m^3 / m^3.}

\item{pk}{air pressure (kPa)}

\item{PAI}{Vector of Plant Area Indices for each canopy layer}

\item{TT}{Cumulative conductivity time to each canopy node (s)}
}
\value{
A list of average values for merged canopy layers of the following components:

\code{tc} air temperature (deg C)

\code{gha} molar heat conductance between leaf and air (mol / m^2 / sec)

\code{gt} molar heat conductance between air layers due to turbulent convection (mol / m^2 / sec)

\code{zla} mean leaf-air distance (m)

\code{z} height above ground of canopy nodes (m)

\code{ph} molar density of air layers (mol / m3)

\code{cp} specific heat of air layer at constant pressure (J / mol / K)

\code{Vo} mole fraction of water vapour in air (mol / mol)

\code{ea} vapur pressure in current time step (kPa)

\code{lambda} Latent heat of vapourization of water (J / mol)

\code{X} Heat to add during time step (deg C)

\code{vden} Volumetric density of vegetation (m^3 / m^3.

\code{m} number of merged canopy layers

\code{PAI} Plant Area Index (m^2 / m^2)

\code{TT} Cumulative conductivity time to each canopy node (s)
}
\description{
Merges canopy airlayers and calculates average molar conductances,
vapour, temperature, leaf area etc
}
