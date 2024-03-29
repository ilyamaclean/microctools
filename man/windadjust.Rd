% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplemodel.R
\name{windadjust}
\alias{windadjust}
\title{Apply height adjustment to wind speed}
\usage{
windadjust(ui, zi, zo)
}
\arguments{
\item{zi}{height of input wind speed (m)}

\item{zo}{height for which wind speed is required (m)}

\item{u}{wind speed at height \code{zi} (m/s)}
}
\value{
wind speed at height \code{zo} (m/s)
}
\description{
Function to apply height adjustment to wind speed,
assuming a reference short grass surface
}
\examples{
zo<-c(20:100)/10
uo<-windadjust(1,2,zo)
plot(zo~uo, type="l", xlab="Wind speed (m/s)", ylab="Height (m)")
}
