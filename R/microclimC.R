#' Calculates radiation absorbed by leaf
#'
#' @description Calculates the flux density of radiation absorbed by leaf.
#'
#' @param Rsw Incoming shortwave radiation (W / m2)
#' @param tme POSIXlt object of time(s)
#' @param tair air temperature (deg C)
#' @param tground ground temperature (deg C)
#' @param lat latitude (decimal degrees)
#' @param long longitude (decimal degrees)
#' @param PAIc Vector of cumulative plant area indices for each canopy layer
#' @param pLAI Proportion of plant area that is green vegetation
#' @param x the ratio of vertical to horizontal projections of leaf foliage
#' @param refls reflectivity of green vegetation to shortwave radiation
#' @param refw reflectivity of woody vegetation to shortwave radiation
#' @param vegem emissivity of vegetation
#' @param skyem sky emissivity
#' @param dp proportion of `Rsw` that is diffuse radiation. If not provided, then calculated using [difprop()]
#' @param merid optional longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst optional value representing the time difference from the timezone meridian
#'
#' @return a list of with the following components:
#' @return `aRsw` absorbed Shortwave radiation (W / m2)
#' @return `aRlw` absorbed longwave radiation (W / m2)
#' @return `ref` mean area-wighted reflectivity of vegetation (green and woody)
#' @export
leafabs <-function(Rsw, tme, tair, tground, lat, long, PAIc, pLAI, x, refls, refw, vegem, skyem, dp = NA,
                 merid = round(long/15, 0) * 15, dst = 0) {
  jd<-julday(tme = tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  if (is.na(dp)) dp<-difprop(Rsw,jd,lt,lat,long,merid=merid,dst=dst)
  sa<-solalt(lt,lat,long,jd,merid=merid,dst=dst)
  ref<-pLAI*refls+(1-pLAI)*refw
  aRsw <- (1-ref) * cansw(Rsw,dp,jd,lt,lat,long,PAIc,x,ref,merid=merid,dst=dst)
  aRlw <- canlw(tair, PAIc, 1-vegem, skyem = skyem)$lwabs
  return(list(aRsw=aRsw, aRlw=aRlw, ref=ref))
}
#' Calculates radiation emitted by leaf
#'
#' @description Calculates the flux density of radiation emitted by the leaf.
#'
#' @param tc temperature (deg C)
#' @param vegem emissivity of leaf
#'
#' @return Flux density of radiation emitted by leaf (W / m2)
#' @export
leafem <- function(tc, vegem) {
 eRlw <- vegem * 5.67*10^-8 * (tc + 273.15)^4
 eRlw
}
#' Thomas algorithm for solving simultanious heat fluxes
#'
#' @description `Thomas` implements the Thomas algorithm for solving simultanious heat
#' fluxes between soil / air layers.
#'
#' @param tc vector of soil and air temperatures (deg C) from previous timestep (see details)
#' @param tmsoil temperature (deg C) of deepest soil layer. Typically mean annual temperature
#' @param tair air temperature at reference height 2 m above canopy in current time step (deg C)
#' @param k vector of thermal conductances between layers (W / m^2 / K) (see details)
#' @param cd thermal heat capacity of layers (W / m^2 / K)
#' @param f forward / backward weighting of algorithm (see details)
#' @param X vector of temperatures to be added resulting from e.g. leaf heat fluxes or radiation
#' absorbed by top soil layer
#' @return a vector of temperatures (deg C) for each soil / air layer for current time step. The first value
#' is `tair` and the last `tmsoil`
#' @export
#' @details The vector `tc` must be ordered with reference air temperature first and the soil temperature
#' of the  deepest layer last. I.e. the length of the vector `tc` is the number of nodes + 2.
#' The vector `k` is the conductivity between each node and that diectly below it, the first value
#' representing conductivity between reference height and the top canopy node. I.e. the length
#' of the vector `k` is the number of nodes + 1. The vector `cd` is the heat storage at each node.
#' I.e. the length of the vector `cd` is the same as the number of nodes. The  weighting factor `f`  may
#' range from 0 to 1. If `f` = 0, the flux is determined by the temperature difference at the beginning
#' of the time step. If `f` = 0.5, the average of the old and new temperatures is used to compute heat flux.
#' If `f` = 1, fluxes are computed using only the new temperatures. The best value to use
#' for `f` is determined by considerations of numerical stability and accuracy and experimentation
#' may be required. If `f` = 0  more heat transfer between nodes is predicted than would actually
#' occur, and can therefore become unstable if time steps are too large. When `f` > 0.5,
#' stable solutions will always be obtained, but heat flux will be underestimated. The
#' best accuracy is obtained with `f` around 0.4, while best stability is at `f` = 1.
#' A typical compromise is `f` = 0.6.
Thomas <- function(tc, tmsoil, tair, k, cd, f = 0.6, X = 0) {
  m <- length(tc) - 2
  tn<-rep(0,m+2)
  tn[m+2]<-tmsoil
  tn[1]<-tair
  g<-1-f
  a <- c(0,0); b <- 0; cc <-0; d <- 0
  xx<-(2:(m+1))
  cc[xx]<- -k[xx]*f
  a[xx+1]<-cc[xx]
  b[xx]<-f*(k[xx]+k[xx-1])+cd
  d[xx]<-X+g*k[xx-1]*tc[xx-1]+(cd-g*(k[xx]+k[xx-1]))*tc[xx]+g*k[xx]*tc[xx+1]
  d[2]<-d[2]+k[1]*tn[1]*f
  d[m+1] <- d[m+1] + k[m+1] * f * tn[m+2]
  for (i in 2:m) {
    cc[i]<-cc[i]/b[i]
    d[i]<-d[i]/b[i]
    b[i+1]<-b[i+1]-a[i+1]*cc[i]
    d[i+1]<-d[i+1]-a[i+1]*d[i]
  }
  tn[m+1] <- d[m+1] / b[m+1]
  for (i in m:2) {
    tn[i]<-d[i]-cc[i]*tn[i+1]
  }
  tn
}
#' Thomas algorithm for solving simultanious vapour fluxes
#'
#' @description `ThomasV` implements the Thomas algorithm for solving simultanious vapour
#' fluxes between air layers.
#'
#' @param Vo a vector of air vapour concentrations for each canopy node in the previos timestep (mol fraction)
#' @param tn vector of air temperatures (deg C) for each canopy node in the current timestep (deg C)
#' @param pk atmospheric pressure (kPa)
#' @param theta Volumetric water content of the upper most soil layer in the current time step (m3 / m3)
#' @param thetap Volumetric water content of the upper most soil layer in the previous time step (m3 / m3)
#' @param relhum relative humidity (percentage) at reference height 2 m above canopy in current time step (percentage)
#' @param tair air temperature at reference height 2 m above canopy in current time step (deg C)
#' @param tsoil temperature of upper soil layer in current time step (deg C)
#' @param zth heightdifference between each canopy node and that directly below it. the first value is
#' the height difference between the lowest canopy node and the ground
#' @param gt vector of molar conductances between each canopy node at that directly below it (mol / m^2 / sec).
#' The first value is the conductivity between the ground and the lowest node, and the last value the
#' conductivity between the highest node and reference height.
#' @param Vflux Total vapour flux from leaves to air (mol /m^3)
#' @param f forward / backward weighting of algorithm (as for [Thomas()])
#' @param previn a list of model outputs form the previous timestep
#' @param soilp a list of soil parameters as returned by [soilinit()]
#' @return a vector of vapour concentrations expressed as mole fractions for each canopy node in the
#' current time step. The first value is that for the ground and the last value that at reference height
#' @export
#' @seealso [Thomas()]
ThomasV <- function(Vo, tn, pk, theta, thetap, relhum, tair, tsoil, zth, gt, Vflux, f = 0.6, previn, soilp) {
  m<-length(zth)
  ph<-phair(tn,pk)
  ea<-0.6108*exp(17.27*tair/(tair+237.3))*(relhum/100)
  eap<-0.6108*exp(17.27*previn$tair/(previn$tair+237.3))*(previn$relhum/100)
  Vair<-ea/pk
  rhsoil<-soilrh(theta,soilp$b,soilp$psi_e,soilp$Smax,tsoil)
  rhsoilp<-soilrh(thetap,soilp$b,soilp$psi_e,soilp$Smax,previn$soiltc[1])
  rhsoil[rhsoil>1]<-1
  rhsoilp[rhsoilp>1]<-1
  eas<-0.6108*exp(17.27*tsoil/(tsoil+237.3))*rhsoil
  easp<-0.6108*exp(17.27*previn$tsoil/(previn$tsoil+237.3))*rhsoilp
  Vsoil<-eas/pk
  Vo<-c(easp/previn$pk,Vo,eap/previn$pk)
  Vn<-Thomas(rev(Vo), Vsoil, Vair, rev(gt), rev(ph), f, Vflux)
  Vn<-rev(Vn)
}
#' Calculates wind profile for entire canopy
#'
#' @description calculates wind speed at any given point above or below canopy
#' from wind speed at reference height
#'
#' @param ui wind at reference height (m / s)
#' @param zi height of wind measurement (m)
#' @param zo height for which wind is required (m)
#' @param a attenuation coefficient as returned by [attencoef()]
#' @param PAI plant area index (m / m)
#' @param hgt height of canopy (m)
#' @param psi_m diabatic correction factor as return by [diabatic_cor()]
#' @param hgtg height of ground vegetation layer below canopy (m)
#' @param zm0 roughness length (m) of vegetation layer below canopy as returned by [roughlength()]
#' @return wind speed (m /s) at height `zo`
#' @export
#' @examples
#' zo<- c(0:200) / 100
#' uz <- windprofile(2, 2, zo, 2.3, 3, 2)
#' plot(zo ~ uz, type = "l", xlab = "wind speed", ylab = "height")
windprofile <- function(ui, zi, zo, a, PAI, hgt, psi_m = 0, hgtg = 0.05 * hgt, zm0 = 0.004) {
  if (length(as.vector(hgt)) == 1) hgt <- hgt + zo * 0
  if (length(as.vector(psi_m)) == 1) psi_m <- psi_m + zo * 0
  if (length(as.vector(a)) == 1) a <- a + zo * 0
  if (length(as.vector(hgtg)) == 1) hgtg <- hgtg + zo * 0
  d <- zeroplanedis(hgt, PAI)
  zm <- roughlength(hgt, PAI, zm0)
  dg <- 0.65 * hgtg
  zmg <- 0.1 * hgtg
  ln1 <- log((zi - d) / zm) + psi_m
  uf <- 0.4 * (ui / ln1)
  # zo above canopy
  ln2 <- suppressWarnings(log((zo - d) / zm) + psi_m)
  uo <- (uf / 0.4) * ln2
  # zo below canopy
  sel <- which(zo < hgt)
  ln2 <- log((hgt[sel] - d[sel]) / zm[sel]) + psi_m[sel]
  uo[sel] <- (uf[sel] / 0.4) * ln2
  # zo above 10% of canopy hgt, but below top of canopy
  sel <- which(zo > (0.1 * hgt) & zo < hgt)
  uo[sel] <- uo[sel] * exp(a[sel] * ((zo[sel] / hgt[sel]) - 1))
  # zo below 10% of canopy hgt
  sel <- which(zo <= (0.1 * hgt))
  uo[sel] <- uo[sel] * exp(a[sel] * (((0.1 * hgt[sel]) / hgt[sel]) - 1))
  ln1 <- log((0.1 * hgt[sel]) / zmg[sel]) + psi_m[sel]
  uf <- 0.4 * (uo[sel] / ln1)
  ln2 <- log(zo[sel] / zmg[sel]) + psi_m[sel]
  uo[sel] <- (uf / 0.4) * ln2
  uo
}
#' Calculates wind profile for individual canopy layers
#'
#' @description Used for calculating wind speed profile in canopies with variable
#' `PAI`
#' @param uh wind speed at top of canopy layer (m /s)
#' @param z height of canopy layer node (m)
#' @param hgt height to top of canopy layer (m)
#' @param PAI Plant Area Index of canopy layer (m / m)
#' @param x the ratio of vertical to horizontal projections of leaf foliage
#' @param lw mean leaf width (m)
#' @param cd drag coefficient
#' @param iw turbulence intensity
#' @param phi_m diabatic correction factor
#' @param edgedist optional numeric value indicating distance (m) to open
#' habitat (see details)
#' @param uref wind at reference height (m) (see details)
#' @param zref height (m) of `uref`
#' @details if 'edgedist' not NA, then horizontal wind component is also added. Here,
#' the wind profile of a reference grass surface (height = 0.12, Plant Area Index = 1.146)
#' is calaculate, and at any given height `z` an attenuated wind speed inside thew canopy
#' calculated, with the degree of attenuation determined by foliage density and distance
#' from the edge. If the the attenuated horizontal wind speed exceeds the wind speed due
#' to vertical attenuation, the horizontal wind speed is returned. The parameter `uref`
#' represents the wind speed at height `zref` above the reference grass surface. If `edgedist`
#' is set to NA (the default) the horizontal wind component and hence `uref` and `zref` are
#' ignored.
#' @return wind speed at height of canopy node (m / s)
#' @export
#' @examples
#' # ==== Generate plant area index values
#' m <- 100
#' hgt <- 10
#' z<-c(1:m) * (hgt / m)
#' PAI <- PAIgeometry(m, 3, 7, 70)
#' plot(z~PAI, type = "l")
#' cPAI <- cumsum(PAI)
#' # ==== Calculate at top of canopy
#' uref <- 2
#' a <- attencoef(hgt, 3, 1)
#' uh <- windprofile(uref, hgt + 2, hgt, a, 3, hgt)
#' # === Calculate canopy profile (near edge)
#' uz1 <- 0
#' for (i in m:1) {
#'   uz1[i] <- windcanopy(uh, z[i], z[i] + 0.05, cPAI[i], edgedist = 5, uref = uref)
#'   uh <- windcanopy(uh, z[i] - 0.05, z[i] + 0.05, cPAI[i], edgedist = 5, uref = uref)
#' }
#' # === Calculate canopy profile (far from edge)
#' uh <- windprofile(uref, hgt + 2, hgt, a, 3, hgt)
#' uz2 <- 0
#' for (i in m:1) {
#'   uz2[i] <- windcanopy(uh, z[i], z[i] + 0.05, cPAI[i], edgedist = 500, uref = uref)
#'  uh <- windcanopy(uh, z[i] - 0.05, z[i] + 0.05, cPAI[i], edgedist = 500, uref = uref)
#' }
#' plot(z ~ uz1, type = "l", xlab = "wind speed", ylab = "height", xlim = c(0,1.5),
#'      col = rgb(1,0,0,0.5), lwd = 2)
#' par(new = TRUE)
#' plot(z ~ uz2, type = "l", xlab = "", ylab = "", xlim = c(0,1.5),
#'      col = rgb(0,0,1,0.5), lwd = 2)
windcanopy <- function(uh, z, hgt, PAI = 3, x = 1, lw = 0.05, cd = 0.2,
                       iw = 0.5, phi_m  = 1, edgedist = NA, uref, zref = hgt + 2) {
  a <- attencoef(hgt, PAI, x, lw, cd, iw, phi_m)
  uz <- uh * exp(a * (z / hgt - 1))
  # horizontal wind component
  if (is.na(edgedist) == F) {
    ah <- attencoef(0.12, 1.146, 0.1, 0.02)
    uhr <- windprofile(uref, zref, z, ah, 1.146, 0.12)
    a2 <- attencoef(hgt, PAI, 1/x, lw, cd, iw, phi_m)
    uhr <- uhr * exp(a2 * ((hgt - edgedist) / hgt - 1))
    uz <- pmax(uz,uhr)
  }
  uz
}
#' Calculates leaf temperature of canopy layers
#'
#' @description Calculates leaf temperature from radiation fluxes and reference air
#' and ground tempertaure
#'
#' @param tair air temperature at two metres above canopy (deg C)
#' @param relhum relative humidity at two metres above canopy (Percentage)
#' @param pk air pressure at two metres above canopy (kPa)
#' @param timestep duration of time step (s)
#' @param z height of canopy layer nodes (m)
#' @param gt heat conductance by turbulent convection (mol / m^2 / s)
#' @param gha leaf-air heat conductance (mol / m^2 / s)
#' @param gv leaf-air vapour conductance (mol / m^2 / s)
#' @param Rabs Flux density of absorbed radiation as returned by [leafabs()] (W/m2)
#' @param previn list of values from previous timestep
#' @param vegp list of vegetation paramaters (see e.g. `vegparams` dataset)
#' @param soilp list Soil paramaters (see e.g. `soilparams` dataset)
#' @param theta volumetric soil moisture fraction of top soil layer (m^3 / m^3)
#'
#' @return a list with the following elements:
#' @return `tn` air temperature of each layer (deg C)
#' @return `tleaf` leaf temperature of each layer (deg C)
#' @return `ea` vapour pressure (kPa)
#' @return `gtt` heat conductance to reference height two m above canopy
#' @return `Rem` emitted radiation (W / m^2)
#' @return `H` Sensible heat flux from leaf to air (W / m^2)
#' @return `L` Latent heat of vapourisation from leaf to air (W / m^2)
#' @return `Lc` Latent heat of condenstation from leaf to air (W / m^2)
#' @return `Vflux` Vapour flux to air from leaves (mol / m^2 / s)
#' @export
#' @details `leaftemp` computes the average leaf and air temperature of each canopy layer based on
#' radiation and evaorative fluxes and reference air and ground temperature. The function
#' automatically determines whether heat storage should be considered, based the specific heat
#' capacity of the vegetation layer and the time step of the model.
#' @examples
#' # Generate paramaters for function:
#' tme <- as.POSIXlt(0, origin = "2020-05-04 12:00", tz = "GMT")
#' previn <- paraminit(20, 10, 10, 15, 80, 11, 500)
#' vegp <- habitatvars(3, 50, -5, tme, m = 20)
#' soilp <- soilinit("Loam")
#' z<-c((1:20) - 0.5) / 20 * vegp$hgt
#' # run function (setting conductances in current time step to same as in previous):
#' ltemp <- leaftemp(11, 80, 101.3, 60, z, previn$gt, previn$gha, previn$gv, previn$Rabs, previn,
#'                   vegp, soilp, 0.3)
#' plot(z ~ ltemp$tleaf, type = "l", xlab = "Leaf temperature", ylab = "Height")
leaftemp <- function(tair, relhum, pk, timestep, z, gt, gha, gv, Rabs, previn, vegp, soilp, theta) {
  lambda <- -42.575*tair+44994
  # Sort out thicknesses
  m<-length(gt)-1
  zth<-c(z[2:m]-z[1:(m-1)],vegp$hgt-(z[m]+z[m-1])/2)
  zla<-mixinglength(vegp$hgt,vegp$PAI,vegp$x,vegp$lw)
  # Sort out conductivitites
  gt<-0.5*gt+0.5*previn$gt
  gv<-0.5*gv+0.5*previn$gv
  gha<-0.5*gha+0.5*previn$gha
  mtref<-0.5*tair+0.5*previn$tair
  mrh<-0.5*relhum+0.5*previn$relhum
  igtt<-rep(1/gt[m+1],m+1)
  igtt2<-1/gt[1]
  for (i in m:1) igtt[i]<-igtt[i+1]+1/gt[i]
  for(i in 2:m) igtt2[i]<-igtt2[i-1]+1/gt[i]
  gtt<-1/igtt[-1]
  gtt2<-1/igtt2
  zref<-(vegp$hgt+2)-z
  mpk<-0.5*previn$pk+0.5*pk
  # Vapour pressure
  esj<-0.6108*exp(17.27*previn$tc/(previn$tc+237.3))
  eaj<-(previn$rh/100)*esj
  estl<-0.6108*exp(17.27*previn$tleaf/(previn$tleaf+237.3))
  esref<-0.6108*exp(17.27*mtref/(mtref+237.3))
  eref<-(mrh/100)*esref
  delta <- 4098*(0.6108*exp(17.27*previn$tleaf/(previn$tleaf+237.3)))/(previn$tleaf+237.3)^2
  rhsoil<-soilrh(theta,soilp$b,-soilp$psi_e,soilp$Smax, previn$soiltc[1])
  esoil<-rhsoil*0.6108*exp(17.27*previn$soiltc[1]/(previn$soiltc[1]+237.3))
  # Test whether steady state
  test<-pmax(timestep*gtt/zref,timestep*gv/zla,timestep*gtt2/z)
  sel<-which(test>1)
  btm<-(1/timestep)+0.5*(gtt/zref+gtt2/z+gv/zla)
  ae<-eaj+0.5*((gtt/zref)*(eref-eaj)+(gtt2/z)*(esoil-eaj)+(gv/zla)*(estl-eaj))/btm
  be<-(0.25*gv*delta)/btm
  tp <- eaj+0.5*((gtt/zref)*eref+(gtt2/z)*esoil+(gv/zla)*estl)
  tp2 <-0.5*(gv/zla)*delta
  btm<-1+0.5*((gtt/zref)+(gtt2/z)+(gv/zla))
  ae[sel]<-tp[sel]/btm[sel]
  be[sel]<-tp2[sel]
  PAIm<-vegp$PAI/zth
  # Air temperature
  # Test whether steady state
  test<-pmax(timestep*gtt/zref,timestep*gha/zla,timestep*gtt2/z)
  # ~Transient
  sel<-which(test>1)
  ph<-phair(previn$tc,previn$pk)
  cp<-cpair(previn$tc)
  vden<-vegp$thickw*vegp$PAI
  ma<-(timestep*PAIm*(1-vden))/cp*ph
  K1<-gtt*cp/zref; K2<-gtt2*cp/z; K3<-gha*cp/zla;
  K4<-(lambda*gv)/(zla*mpk); K5<-(lambda*gtt2)/(z*mpk)
  btm<-1+0.5*ma*(K1+K2+K3)
  aL<-(previn$tc+0.5*ma*(K1*mtref+K2*previn$soiltc[1]+K3*previn$tleaf+K4*(estl-ae)+K5*(esoil-ae)))/btm
  bL<-ma*(0.25*K3+0.25*K4*delta-0.5*be*(K4+K5))/btm
  # ~Steady state
  K1<-gtt[sel]*cp[sel]; K2<-gtt2[sel]*cp[sel]; K3<-gha[sel]*cp[sel]
  K4<-lambda*gv[sel]/mpk; K5<-lambda*gtt2[sel]/mpk
  aL[sel]<-(K1*mtref+K2*previn$soiltc[1]+K3*previn$tleaf[sel]+K4*(estl[sel]-ae[sel])+K5*(esoil-ae[sel]))/(K1+K2+K3)
  bL<-(0.5*(K3+K4*delta[sel]-be[sel]*(K4+K5)))/(K1+K2+K3)
  # Sensible Heat
  aH<-cp*gha*(previn$tleaf-aL)
  bH<-cp*gha*(0.5-0.5*bL)
  # Latent Heat
  aX<-(lambda*gv/mpk)*(estl-ae)
  bX<-(lambda*gv/mpk)*(delta-be)
  # Radiation
  Rabs<-0.5*Rabs+0.5*previn$Rabs
  sb<-5.67*10^-8
  aR<-vegp$vegem*sb*(previn$tleaf+273.15)^4
  bR<-vegp$vegem*sb*2*(previn$tleaf+273.15)^3
  # Leaf temperature
  Ch<-vden*vegp$cpw*vegp$phw
  ml<-timestep*PAIm/Ch
  dTL<-(ml*(Rabs-aR-aX-aH))/(zla+ml*(bR+bX+bH))
  # Latent heat if condensation
  tn<-aL+bL*dTL
  ea<-ae+be*dTL
  es<-0.6108*exp(17.27*tn/(tn+237.3))
  tn2<-(237.3*log(es/0.6108))/(17.27-log(es/0.6108))
  sel<-which(ea>es)
  Lc<-lambda*((ea-es)/pk)*ph
  dTL2<-(ml/zla)*Lc
  # Temperatures and fluxes
  tn[sel]<-tn2[sel]
  dTL[sel]<-dTL[sel]+dTL2[sel]
  Lc2<-ea*0; Lc2[sel]<-Lc[sel]
  eam<-ae+be*dTL
  neaj<-eaj+2*(eam-eaj)
  # Vapour flux from leaves
  mTl<-previn$tleaf+0.5*dTL
  mVl<-(0.6108*exp(17.27* mTl/(mTl+237.3)))/mpk
  mVa<-eam/mpk
  Vflux<-(gv/zla)*(mVl-mVa)
  # Save outputs
  return(list(tn=tn, tleaf=previn$tleaf+dTL, ea=neaj, gtt=gtt,
              Rem=aR+bR*dTL, H=aH+bH*dTL, L=aX+bX*dTL, Lc=Lc2, Vflux = Vflux))
}
#' Initialise paramaters for first time step of model
#'
#' @description Generates a set of climate and conductivity parameters for running the
#' first time step of the model
#'
#' @param m number of canopy layer nodes
#' @param sm number of soil layers
#' @param hgt height of canopy (m)
#' @param tair air temperature at 2 m above canopy (deg C)
#' @param relhum relative humidity at 2 m above canopy (percentage)
#' @param tsoil temperature of deepest soil layer. Usually ~mean annual temperature
#' (deg C). See details.
#' @param Rsw total incoming shortwave radiation (W / m^2)
#' @return a list with the following elements:
#' @return `tc` a vector of air temperatures for each canopy layer (deg C)
#' @return `soiltc` a vector of airsoil temperatures for each soil layer (deg C)
#' @return `tleaf` a vector of leaf temperatures for each canopy layer (deg C)
#' @return `rh` a vector of relative humidities
#' @return `relhum` relative humidity at 2 m above canopy (percentage)
#' @return `tair` air temperature at 2 m above canopy (deg C)
#' @return `pk` pressure at 2 m above canopu (kPA)
#' @return `Rabs` Absorbed radiation (W / m^2)
#' @return `gt` Conductivity in air of each canopy layer node (mol/m^2/sec)
#' @return `gv` Leaf conductivity to vapour loss for each canopy layer node  (mol/m^2/sec)
#' @return `gha` Conductivity between air and leaf for each canopy layer node (mol/m^2/sec)
#'
#' @importFrom stats spline
#' @export
#' @examples
#' paraminit(20, 10, 10, 15, 80, 11, 500)
#'
#' @details All values are approximate. Values for `tc` and `tsoil` are derived by
#' linear intepolation between `tair` and `tsoil`. Values for `Rabs` are derived from
#' `Rsw` but attenuate through the canopy. Values for `tleaf` are derived
#' from `tc` and `Rabs`. Values for `rh` are the same as `relhum`. Values for `gt`
#' `gv` and `gha` are typical for decidious woodland with wind above canopy at 2 m/s.
#' `gt` is scaled by canopy height and `m` (and hence distance between nodes). The first
#' value represents conductivity between the ground and the lowest canopy node. The last
#' value represents conductivity between the air at 2 m above canopy and the highest
#' canopy node.
paraminit <- function(m, sm, hgt, tair, relhum, tsoil, Rsw) {
  tcb <- spline(c(1,2), c(tsoil, tair), n = m + sm)
  tc <- tcb$y[(sm + 1): (m + sm)]
  soiltc <- rev(tcb$y[1:sm])
  rh <- rep(relhum, m)
  Rabs <- spline(c(1, 2), c(0.2 * Rsw + 20, Rsw + 20), n = m)$y
  tleaf <- tc + 0.01 * Rabs
  gt <- rep(50, m + 1) * (m / hgt) * (2 / m)
  gt[1] <- 2 * gt[1]
  gt[m+1] <- gt[m+1] * (hgt / m) * 0.5
  gv <- spline(c(1, 2), c(0.25, 0.32), n = m)$y
  gha <- spline(c(1, 2), c(0.13, 0.19), n = m)$y
  z<-c((1:m)-0.5)/m*hgt
  sz<-2/sm^2.42*c(1:m)^2.42
  return(list(tc = tc, soiltc = soiltc, tleaf = tleaf, z = z, sz = sz, rh = rh, relhum = relhum,
              tair = tair, tsoil = tsoil, pk = 101.3, Rabs = Rabs, gt = gt, gv = gv, gha = gha))
}
#' Returns soil parameters for a given soil type
#'
#' @description Returns soil parameters needed to run the microclimate model
#' for a given soil type
#'
#' @param soiltype one of `Sand`, `Loamy sand`, `Sandy loam`, `Loam`, `Silt`,
#' `Silt loam`, `Sandy clay loam`, `Clay loam`, `Silty clay loam`, `Sandy clay`,
#' `Silty clay` or `Clay`.
#'
#' @return a list with the following items:
#' @return `Soil.type` description of soil type
#' @return `Smax` Volumetric water content at saturation (m^3 / m^3)
#' @return `Smin` Residual water content (m^3 / m^3)
#' @return `alpha` Shape parameter of the van Genuchten model (cm^-1)
#' @return `n` Pore size distribution parameter (dimensionless, > 1)
#' @return `Ksat` Saturated hydraulic conductivity (cm / day)
#' @return `Vq` Volumetric quartz content of soil
#' @return `Vm` Volumetric mineral content of soil
#' @return `Vo` Volumetric organic content of soil
#' @return `Mc` Mass fraction of clay
#' @return `rho` Soil bulk density (Mg / m^3)
#' @return `b` Shape parameter for Campbell model (dimensionless, > 1)
#' @return `psi_e` Matric potential (J / m^3)
#'
#' @export
#' @examples
#' soilinit("Loam")
soilinit <- function(soiltype) {
  sel <- which(soilparams$Soil.type == soiltype)
  soilp <- soilparams[sel,]
  as.list(soilp)
}
#' Run canopy model for a single time step
#' @description run the below-canopy model for a single timestep
#'
#' @param climvars a of climate variables needed to run the run the model for one timestep (see dataset [climvars()])
#' @param previn a list of model outputs form the previous timestep as returned initially by [paraminit()]
#' @param vegp a list of vegetation parameters as returned by [habitatvars()]
#' @param soilp a list of soil parameters as returned by [soilinit()]
#' @param timestep length of model timestep (s)
#' @param tme POSIXlt object of the date and time of the current time step
#' @param lat latitude of location (decimal degrees)
#' @param long lonitude of location (decimal degrees)
#' @param edgedist distance to open ground (m)
#' @param reqhgt optional height for which temperature is required (see details)
#' @param sdepth depth of deepest soil node (m)
#' @param zu height above ground of reference climate measurements (m)
#' @param theta volumetric water content of upper most soil layer in current time step (m^3 / m^3)
#' @param thetap volumetric water content of upper most soil layer in previous time step (m^3 / m^3)
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT).
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if `merid` = 0).
#' @param n forward / backward weighting for Thomas algorithm (see [Thomas()])
#' @return a list of of model outputs for the current timestep with the same format as `previn`
#' @export
#' @details model outputs are returned for each canopy node, with the number of canopy nodes (m)
#' determined by `previn`. Canopy nodes are spaced at equal heights throughout the canopy as in the
#' example. If `reqhgt` is set, the canopy node nearest to that height is set at the value specified.
#' @examples
#' # Create initail parameters
#' tme <- as.POSIXlt(0, origin = "2020-05-04 12:00", tz = "GMT")
#' previn <- paraminit(20, 10, 10, 15, 80, 11, 500)
#' vegp <- habitatvars(4, 50, -5, tme, m = 20)
#' z<-c((1:20)-0.5)/20*vegp$hgt
#' soilp<- soilinit("Loam")
#' climvars <- list(tair=16,relhum=90,pk=101.3,u=2.1,tsoil=11,skyem=0.9,Rsw=500,dp=NA,
#'                  psi_h=0,psi_m=0,phi_m=0)
#' # Run model 100 times for current time step
#' for (i in 1:100) {
#'   plot(z ~ previn$tc, type = "l", xlab = "Temperature", ylab = "Height", main = i)
#'   previn <- runcanopy(climvars, previn, vegp, soilp, 60, tme, 50, -5)
#' }
runcanopy <- function(climvars, previn, vegp, soilp, timestep, tme, lat, long, edgedist = 100,
                      sdepth = 2, reqhgt = NA, zu = 2, theta = 0.3, thetap = 0.3, merid = 0,
                      dst = 0, n = 0.6) {
  # =============   Unpack climate variables ========== #
  m <- length(previn$tc)
  tair<-climvars$tair; relhum<-climvars$relhum; pk<-climvars$pk; u<-climvars$u
  tsoil<-climvars$tsoil; skyem<-climvars$skyem; Rsw<-climvars$Rsw; dp<-climvars$dp
  psi_h<-climvars$psi_h; phi_m<-climvars$phi_m; psi_m<-climvars$psi_m
  # ========== Calculate baseline variables ============== #
  tc<-previn$tc; ppk<-previn$pk; hgt<-vegp$hgt
  ph<-phair(tc,ppk); pha<-phair(tair,pk) # molar density of air
  cp<-cpair(tc) # specific heat of air
  lambda <- -42.575*tc+44994 # Latent heat of vapourisation (J / mol)
  # Adjust wind to 2 m above canopy
  u2<-u*log(67.8*hgt-5.42)/log(67.8*zu-5.42)
  # Generate heights of nodes
  z<-c((1:m)-0.5)/m*vegp$hgt
  if (is.na(reqhgt) == F) z[abs(z-reqhgt)==min(abs(z-reqhgt))][1]<-reqhgt
  zt<-z[2:(m)]-z[1:(m-1)] # difference in height between layers
  # ========== Calculate wind speed and turbulent conductances ======== #
  ac<- attencoef(hgt,sum(vegp$PAI),vegp$x,vegp$lw,vegp$cd,vegp$iw,phi_m)
  uz<-rep(0,m); gt<-rep(0,m)
  uh<-windprofile(u2,hgt+2,hgt,ac[m],sum(vegp$PAI),hgt,psi_m,vegp$hgtg,vegp$zm0)
  uz[m]<-windcanopy(uh,z[m],hgt,sum(vegp$PAI),vegp$x,vegp$lw,vegp$cd,vegp$iw[m],
                    phi_m,edgedist,u,zu)
  gt[m]<-gcanopy(uh,z[m],z[m-1],tc[m],tc[m-1],hgt,sum(vegp$PAI),vegp$x,vegp$lw,
                 vegp$cd,vegp$iw[m],phi_m,pk)
  tc<-c(previn$soiltc[1],tc)
  # ========== Calculate canopy turbulences ========== #
  for (i in (m-1):1) {
    zi<-z[i+1]+0.5*zt[i]; zo<-z[i+1]-0.5*zt[i]
    uh<-windcanopy(uh,zo,zi,sum(vegp$PAI[1:i]),vegp$x,vegp$lw,vegp$cd,vegp$iw[i],
                   phi_m,edgedist,u,zu)
    if (i>1) {
      z0<-z[i-1]
    } else z0<-0
    uz[i]<-windcanopy(uh,zi,zo,sum(vegp$PAI[1:i]),vegp$x,vegp$lw,vegp$cd,
                      vegp$iw[i],phi_m,edgedist,u,zu)
    gt[i]<-gcanopy(uh,z[i],z0,tc[i+1],tc[i],zi,sum(vegp$PAI[1:i]),vegp$x,vegp$lw,
                   vegp$cd,vegp$iw[i],phi_m,pk)
  }
  # Turbulent air conductivity and layer merge
  gt[m+1]<-gturb(u2,hgt+2,hgt+2,hgt,hgt,sum(vegp$PAI),tair,psi_m,psi_h,vegp$zm0,pk)
  lmm<-layermerge(z,gt,hgt,timestep,ph)
  mrge<-lmm$mrge; gtx<-lmm$gtx; zth<-lmm$zth
  tc <- tc[-1]
  # ========== Calculate absorbed radiation =========== #
  PAIc<-rev(cumsum(rev(vegp$PAI)))
  Rabss<-leafabs(Rsw,tme,tair,previn$soiltc[1],lat,long,PAIc,vegp$pLAI,vegp$x,vegp$refls,
                 vegp$refw,vegp$vegem,skyem,dp,merid,dst)
  Rabs<-Rabss$aRsw+Rabss$aRlw
  # ============= Conductivities =============== #
  # Vapour conductivity
  gv<-layercond(Rabss$aRsw/(1-Rabss$ref),vegp$gsmax,vegp$q50)
  # Leaf conductivity
  tc2<-c(tc[2:length(tc)],tair); dtc<-tc2-tc
  gha<-1.41*gforcedfree(vegp$lw*0.71,uz,tc,dtc,pk)
  tln<-leaftemp(tair,relhum,pk,timestep,z,gt,gha,gv,Rabs,previn,vegp,soilp,theta)
  tleafm <- 0.5*tln$tleaf+0.5*previn$tleaf
  eaj<-0.6108*exp(17.27*tc/(tc+237.3))*(previn$rh/100)
  Vo<-eaj/previn$pk
  # =============== Soil conductivity =========== #
  sm<-length(previn$soiltc)
  cdk<-soilk(timestep,sm,sdepth,theta,soilp$Vm,soilp$Vq,soilp$Mc,soilp$rho)
  sz<-cdk$z
  # conductivity and specific heat
  vden<-(vegp$PAI*vegp$thickw)
  mult<-1-vden
  zla <- mixinglength(vegp$PAI,vegp$hgt,vegp$x,vegp$lw)
  L<-tln$L-tln$Lc
  TT<-cumsum((ph/gt[1:m])*(z-c(0,z[1:(m-1)])))
  lav<-layeraverage(lmm,tc,tleafm,vegp$hgt,gha,gt,zla,z,Vo,tln$Vflux,L,tln$H,vden,ppk,vegp$PAI,TT)
  cda<-lav$cp*lav$ph*(1-lav$vden)*(lav$z[3:(lav$m+2)]-lav$z[1:(lav$m)])/2*timestep
  ka<-lav$gt*c(lav$cp,cpair(previn$tair))
  ka[1:lav$m]<-ifelse(ka[1:lav$m]>cda,cda,ka[1:lav$m])
  k<-c(rev(ka),cdk$k)
  cd<-c(rev(cda),cdk$cd)
  # Heat to add
  ma<-(timestep*lav$PAI)/(c(lav$ph,phair(previn$tair,previn$pk))*
                            c(lav$cp,cpair(previn$tair))*c((1-lav$vden),1))
  Xa<-rev(c(ma*(lav$H+lav$L))); Xa<-Xa[-1]
  ref<-vegp$pLAI[1]*vegp$refls+(1-vegp$pLAI[1])*vegp$refw
  Xs<-(timestep/cdk$cd[1])*(1-vegp$refg)*Rabss$aRsw[1]/(1-ref)
  X<-c(Xa,Xs,rep(0,(sm-1)))
  # Heat exchange between layers
  tc2<-c(previn$tair,rev(lav$tc),previn$soiltc, previn$tsoil)
  tn2<-Thomas(tc2, tsoil, tair, k, cd, n, X)
  tnair<-rev(tn2[1:(lav$m+2)])
  tnsoil<-tn2[(lav$m+2):(length(tn2)-1)]
  # vapour exchange
  zth2<-lav$z[2:(lav$m+1)]-lav$z[1:lav$m]
  Vmflux<-(timestep*lav$Vflux*lav$PAI[1:lav$m])/zth2
  gtmx<-c(lav$ph/zth2*(2*timestep),Inf)
  gt2<-ifelse(lav$gt>gtmx,gtmx,lav$gt)
  tn2<-tnair[-1]; tn2<-tn2[-length(tn2)]
  Vn<-ThomasV(lav$Vo,tn2,pk,theta,thetap,relhum,tair,tnsoil[1],zth2,gt2,Vmflux,n,previn,soilp)
  # Interpolate
  TX<-TT[length(TT)]+(pha/gt[m+1])*2
  T2<-c(0,lav$TT,TX)
  tn<-layerinterp(T2, TT, tnair)
  vn<-layerinterp(T2, TT, Vn)
  # relative humidity
  ea<-vn*pk
  es<-0.6108*exp(17.27*tn/(tn+237.3))
  rh<-(ea/es)*100
  rh[rh>100]<-100
  dataout<-list(tc=tn,soiltc=tnsoil,tleaf=tln$tleaf,z=z,sz=sz,rh=rh,
                relhum=relhum,tair=tair,tsoil=tsoil,pk=pk,Rabs=Rabs,
                gt=gt,gv=gv,gha=gha)
  return(dataout)
}
