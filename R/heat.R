#' Calculate molar density of air
#' @param tc temperature (deg C)
#' @param pk atmospheric pressure (kPa)
#' @return molar density of air (mol / m3)
#' @export
phair <- function(tc = 15, pk = 101.3) {
  tk <- tc + 273.15
  ph <- 44.6 * (pk / 101.3) * (273.15 / tk)
  ph
}
#' Calculate molar specific heat of air
#' @param tc temperature (deg C)
#' @return specific heat of air at constant pressure (J / mol / K)
#' @export
cpair <- function(tc) {
  cp <-  2e-05 * tc^2 + 0.0002 * tc + 29.119
  cp
}
#' Calculate zero plane displacement
#'
#' @param hgt canopy height (m)
#' @param PAI Plant area index
#' @return zero plane displacement (m)
#' @export
zeroplanedis <- function(hgt, PAI = 3) {
  PAI[PAI < 0.1] <- 0.1
  PAI[PAI > 20] <- 20
  m <- 0.0609 * log(PAI) + 0.5894
  d <- m * hgt
  d
}
#' Calculate roughness length governing momentum transfer
#'
#' @param hgt canopy height (m)
#' @param PAI Plant area index
#' @param zm0 minimum roughness length (relfecting ground surface roughness below canopy)
#' @return momentum roughness length (m)
#' @export
roughlength <- function(hgt, PAI = 3, zm0 = 0.004) {
  m <- 0.1 * PAI + 0.08
  sel <- which(PAI > 0.6)
  m[sel] <- -0.0239 * log(PAI[sel]) + 0.1275
  zm <- m * hgt
  zm <- ifelse(zm < zm0, zm0, zm)
  zm
}
#' Calculates forced or free laminer conductance
#'
#' @param d chacteristic dimension of surface (m)
#' @param u wind speed (m/s)
#' @param tc temperature (deg C)
#' @param dtc estimate of temperature differences of surface and air, e.g.
#' from previous time step (see details)
#' @param pk atompsheric pressure (KPa). used for calculating molar density
#' of air.
#' @return conductance (mol / m^2 / sec)
#' @export
#'
#' @details Calculates conductance under forced and free convection and selects
#' whichever is greater. For conductance under free convection an estimate of
#' the temperature difference between the air and  surface is needed (usually
#' from the previous timestep)
#'
#' @examples
#' # As function of `d`
#' d <- c(1:1000) / 100
#' g <- gforcedfree(d, 2, 11, 2)
#' plot(g~d, type = "l")
#' # As function of `u`
#' u <- c(1:200)/100
#' g<-gforcedfree(0.1,u,11,4)
#' plot(g~u, type = "l")
gforcedfree <- function(d, u, tc, dtc, pk = 101.3) {
  # Reynolds etc
  tk <- tc + 273.15
  v <- (0.0908 * tk - 11.531) / 1000000
  Dh <- (0.1285 * tk - 16.247) / 1000000
  Gr <- (9.807  * d^3  * abs(dtc - 0)) / (tk * v^2)
  re <-  (u * d) / v
  Pr <- v / Dh
  # Forced
  ph <- phair(tc, pk)
  cp <-  cpair(tc)
  gh <- (0.34 * Dh * ph * re^0.5 * Pr ^ (1/3)) / d
  # Free
  m <- ifelse(dtc > 0, 0.54, 0.26)
  ghfr <- (m * ph * Dh * (Gr * Pr) ^ (1/4)) / d
  # Set to wichever is higher
  gh <- ifelse(gh > ghfr, gh, ghfr)
  gh
}
#' Calculate mixing length for canopy air transport
#'
#' @param hgt canopy height (m)
#' @param PAI plant area index
#' @param x the ratio of vertical to horizontal projections of leaf foliage
#' @param lw mean leaf width (m)
#' @return mixing length (m). See details
#' @export
#' @details The mixing length is us used to calculate tubulent air transport
#' inside vegetated canopies. It is the mean free pathway, equivelent to the
#' mean leaf-air distance. The leaf distribution angle is used to asses whether
#' leaves are more like grasses or squares (assumed to be like squares if less than
#' 1 and grasses if more than one).
mixinglength <- function(hgt, PAI, x, lw) {
  Ld <- PAI / hgt
  lmg <- (4 * lw / (pi * Ld))^0.5
  lms <- (6 * lw^2 * hgt / (pi * PAI))^(1 / 3)
  wgt<-x/2
  sel <- which(x > 1)
  wgt[sel] <- 1-(0.5/x[sel])
  l_m <- wgt * lms + (1 - wgt) * lmg
  l_m
}
#' Calculate mixing length for canopy air transport
#'
#' @description Calculates mixing length for canopy air transport - equivelent to mean
#' leaf spacing.
#'
#' @param hgt canopy height (m)
#' @param PAI plant area index
#' @param x leaf distribution angle coefficient
#' @param lw mean leaf width (m)
#' @param cd drag coefficient
#' @param iw relative turbulence intensity
#' @param phi_m diabatic correction factor
#' @return dimensionless attenuation coefficient for calculating the wind height profile
#' in vegetation canopies.
#' @export
attencoef <- function(hgt, PAI = 3, x = 0.5, lw = 0.05, cd = 0.2, iw = 0.5, phi_m  = 1) {
  l_m <- mixinglength(hgt, PAI, x, lw)
  a <- ((cd * PAI * hgt) / (2 * l_m * iw))^0.5
  a <- a * phi_m^0.5
  a
}
#' Calculate conductance under turbulent convection above canopy
#'
#' @param u wind speed at height `zu` (m/s)
#' @param zu height of wind speed measurement (m)
#' @param z1 upper height to which conductance is wanted (m)
#' @param z0 lower height from which conductance is wanted (set to heat
#' exchange surface of canopy if NA)
#' @param hgt height of canopy (m)
#' @param PAI plant area index for determining canopy roughness length
#' @param tc temperature used for calculating molar density of air (deg C)
#' @param psi_m diabatic correction factor for momentum transfer
#' @param psi_h diabataic correction factor for heat transfer
#' @param zm0 minimum surface roughness
#' @param pk atmospheric pressure (kPA) used for calculating molar density of air
#' @return molar conductance under turbulent convection above canopy (mol/m^2/s)
#' @export
gturb <- function(u, zu = 2, z1, z0 = NA, hgt, PAI= 3, tc = 15,
                  psi_m = 0, psi_h = 0, zm0 = 0.004, pk = 101.3) {
  d <- zeroplanedis(hgt, PAI)
  zm <- roughlength(hgt, PAI, zm0)
  zh <- 0.2 * zm
  if (is.na(z0)) z0 <- d + zh
  ph <- phair(tc, pk)
  cp <- cpair(tc)
  ustr <- (0.4 * u) / (log((zu - d) / zm) + psi_m)
  g <- (0.4 * ph * ustr) / (log((z1 - d) / (z0 - d)) + psi_h)
  Dh <- (0.1285 * (tc + 273.15) - 16.247) / 1000000
  zref <- d + zh
  gmin <- (ph * Dh) / (z1 - z0)
  g <- ifelse(g < gmin, gmin, g)
  g[is.na(g)] <- gmin
  g
}
#' Calculates conductance under turbulent convection within the canopy
#'
#' @param uh wind speed at height of canopy top (m/s) as returned by [microclimc::windprofile()]
#' @param z1 upper height to which conductance is wanted (m)
#' @param z0 lower height from which conductance is wanted (m)
#' @param tc1 temperature of upper layer (dec C) usually in previous timestep
#' @param tc0 temperature of lower layer (dec C) usually in previous timestep
#' @param hgt height of canopy (m)
#' @param PAI plant area index formixing and attenuation coefficients
#' @param x the ratio of vertical to horizontal projections of leaf foliage
#' @param lw mean leaf width (m)
#' @param cd drag coefficient
#' @param iw relative turbulence intensity
#' @param phi_m diabatic correction factor
#' @param pk atmospheric pressure used for calculating moar density of air (kPA)
#' @return conductance under tubulent convection within the canopy (mol/m^2/sec)
#' @export
gcanopy <- function(uh, z1, z0, tc1, tc0, hgt, PAI = 3, x = 0.5, lw = 0.05,
                    cd = 0.2, iw = 0.5, phi_m = 1, pk = 101.3) {
  a <- attencoef(hgt, PAI, x, lw, cd, iw, phi_m)
  l_m <- mixinglength(hgt, PAI, x, lw)
  tcm <- (tc1 + tc0) / 2
  ph <- phair(tcm, pk)
  tp <- (l_m * iw * ph * uh * a) / phi_m
  e0 <- exp(-a * (z0 / hgt - 1))
  e1 <- exp(-a * (z1 / hgt - 1))
  g <- tp / (e0 - e1)
  # Set minimum
  cp <- cpair(tcm)
  dTdz <- abs((tc1 - tc0) / (z1 - z0))
  Kmin <- (iw * l_m^2 * 3 / 0.74)^ (2/3) *  ((4.6 * dTdz) / (tcm + 273.15))^0.5
  gmin <- (Kmin * ph) / (z1 - z0)
  g <- ifelse(g < gmin, gmin, g)
  g[is.na(g)] <- gmin
  g
}
#' Calculates diabatic correction factor above canopy
#'
#' @description Calculates the diabatic correction factors used in adjustment of
#' wind profiles and calculation of turbulent conductivity above canopy
#'
#' @param tc temperature
#' @param pk atmospheric pressure (kPa)
#' @param H Heat flux (W / m^2)
#' @param uf friction velocity (m/s)
#' @param zi height to which correction factor is wanted (m)
#' @param d zero plane displacement height as returned by [zeroplanedis()]
#' @return a list with the following components:
#' @return `psi_m` diabatic correction factor for momentum transfer
#' @return `psi_h` diabatic correction factor for heat transfer
#'
#' @export
diabatic_cor <- function(tc, pk = 101.3, H = 0, uf, zi = 2, d) {
  Tk <- tc + 273.15
  ph <- phair(tc, pk)
  cp <-  cpair(tc)
  st <- -(0.4 * 9.81 * (zi - d) * H) / (ph * cp * Tk * uf^3)
  # Stable flow
  sel <- which(st < 0) # unstable
  # Stable
  psi_h <- suppressWarnings(6 * log(1 + st))
  psi_m <- psi_h
  # Unstable
  psi_h[sel] <-   -2 * log((1 + (1 - 16 * st[sel])^0.5) / 2)
  psi_m[sel] <- 0.6 * psi_h[sel]
  psi_m <-ifelse(psi_m > 5, 5, psi_m)
  psi_h <-ifelse(psi_h > 5, 5, psi_h)
  return(list(psi_m = psi_m, psi_h = psi_h))
}
#' Calculates diabatic correction factor in canopy
#'
#' @description Calculates the diabatic correction factors used in adjustment of
#' wind profiles and calculation of turbulent conductivity within the canopy
#'
#' @param tc vector of air temperatures at canopy nodes (deg C)
#' @param uz vector of wind speeds at canopy nodes (m/s)
#' @param z vector of heights of canopy nodes (m)
#' @param x leaf angle coefficient
#' @param lw mean leaf width
#' @param PAI vector of plant area index values at each canopy node
#' @return a list with the following components:
#' @return `phi_m` diabatic correction factor for momentum transfer
#' @return `phi_h` diabatic correction factor for heat transfer
#'
#' @export
diabatic_cor_can <- function(tc, uz, z, PAI, x = 1, lw = 0.05) {
  dtc<-tc[2:length(tc)]-tc[1:(length(tc)-1)]
  dz<-z[2:length(z)]-z[1:(length(z)-1)]
  dtdz<-dtc/dz
  tk <- (tc[2:length(tc)]+tc[1:(length(tc)-1)])/2+273.15
  PAIl<-(PAI[2:length(PAI)]+PAI[1:(length(PAI)-1)])/2
  l_m <- mixinglength(dz, PAIl, x, lw)
  u<-(uz[2:length(uz)]+uz[1:(length(uz)-1)])/2
  Ri <- (9.81 / tk) * dtdz * (l_m / u)^2
  Ri[Ri > 0.15] <- 0.15
  st <- (0.74 * (1 + 8.926 * Ri) ^ 0.5 + 2 * 4.7 * Ri - 0.74) / (2 * 4.7 * (1 - 4.7 * Ri))
  sel <- which(dtdz <= 0) # unstable
  st[sel]<-Ri[sel]
  # Stable
  phi_m <- 1 + (6 * st) / (1 + st)
  phi_h <- phi_m
  # Unstable
  phi_m[sel] <- 1 / (1 - 16 * st[sel])^0.25
  phi_h[sel] <- phi_m[sel]^2
  return(list(phi_m = mean(phi_m), phi_h = mean(phi_h)))
}
