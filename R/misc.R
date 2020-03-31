#' Determines whether canopy air layers should be merged
#'
#' @description Determines whether the time step and heat conductance between specified air
#' layers within the canopy are such that temperatures would equalise and hence
#' layers should be merged.
#'
#' @param z vector of heights representing nodes at the mid-point of each canopy layer (m)
#' @param gt conductance by tubulent heat transfer between each node and that below it (see details) (mol / m^2 / s).
#' @param hgt canopy height (m)
#' @param timestep time step of model (s)
#' @param ph molar density of air as returned by [phair()] (mol / m^3)

#' @details
#' The first value of `gt` is conductance between the ground and the lowest node.
#' Subsequent values are for between nodes, as returned by [gcanopy()].
#' The final value is the conductance between the heighest node and a point 2 m above the
#' canopy representing conductance in series between the heighest node and the top of
#' the canopy and the top of the canopy and the air above it, the latter returned
#' by [gturb()].
#'
#' @return a list of with the following components:
#' @return `mrge` a vector of integers for each layer. Those with the same number represent
#' layers that should be merged.
#' @return `zth` thickness of each (unmerged layer)
#' @return `u` vector of unique layer after merging
#' @export
layermerge<-function(z, gt, hgt, timestep, ph = rep(42.24, length(z))){
  m<-length(z)
  zth<-(z[1]+z[2])/2
  for (i in 2:(m-1)) zth[i]<-(z[i]+z[i+1])/2-(z[i]+z[i-1])/2
  zth[m]<-hgt-(z[m]+z[m-1])/2
  zthi<-zth[m]; gtx<-1/gt[m+1]
  mrge<-c(1:(m+1))
  for (i in m:1) {
    gmx<-(zthi*ph[i])/timestep; gtx<-1/((1/gt[i])+(1/gtx))
    zthi<-ifelse(gmx<gtx,zthi+zth[i],zth[i])
    mrge[i]<-ifelse(gmx<gtx,mrge[i+1],mrge[i+1]-1)
    gtx<-ifelse(gmx<gtx,gtx,Inf)
  }
  u<-unique(mrge)
  return(list(mrge=mrge, zth=zth, u=u))
}
#' Merges canopy air layers and calculates average properties
#'
#' @description Merges canopy airlayers and calculates average molar conductances,
#' vapour, temperature, leaf area etc
#'
#' @param lmm list of layers to merge as returned by [layermerge()]
#' @param tc air temperature of unmerged canopy layers (dec C)
#' @param tleaf leaf temperature of unmerged canopy layers
#' @param hgt height of canopy (m)
#' @param gha molar heat conductance between leaf and air (mol / m^2 / sec)
#' @param gt molar heat conductance between air layers due to turbulent convection (mol / m^2 / sec)
#' @param zla mean leaf-air distance (m)
#' @param z height above ground of canopy nodes (m)
#' @param Vo mole fraction of water vapour in air (mol / mol)
#' @param Vflux molar flux density of water vapour from leaves to air (mol / m^2/ sec)
#' @param L Latent heat flux from canopy layer (W / m^2)
#' @param H Sensible heat flux from canopy layer (W / m^2)
#' @param vden  Volumetric density of vegetation (m^3 / m^3.
#' @param pk air pressure (kPa)
#' @param PAI Vector of Plant Area Indices for each canopy layer
#' @param TT Cumulative conductivity time to each canopy node (s)
#' @return A list of average values for merged canopy layers of the following components:
#' @return `tc` air temperature (deg C)
#' @return `tleaf` leaf temperature (dec C)
#' @return `gha` molar heat conductance between leaf and air (mol / m^2 / sec)
#' @return `gt` molar heat conductance between air layers due to turbulent convection (mol / m^2 / sec)
#' @return `zla` mean leaf-air distance (m)
#' @return `z` height above ground of canopy nodes (m)
#' @return `ph` molar density of air layers (mol / m3)
#' @return `cp` specific heat of air layer at constant pressure (J / mol / K)
#' @return `Vo` mole fraction of water vapour in air (mol / mol)
#' @return `Vflux` molar flux density of water vapour from leaves to air (mol / m^2/ sec)
#' @return `lambda` Latent heat of vapourization of water (J / mol)
#' @return `L` Latent heat flux from canopy layer (W / m^2)
#' @return `H` Sensible heat flux from canopy layer (W / m^2)
#' @return `vden` Volumetric density of vegetation (m^3 / m^3.
#' @return `m` number of merged canopy layers
#' @return `PAI` Plant Area Index (m^2 / m^2)
#' @return `TT` Cumulative conductivity time to each canopy node (s)
#' @export
layeraverage<-function(lmm, tc, tleaf, hgt, gha, gt, zla, z, Vo, Vflux, L, H, vden, pk, PAI, TT) {
  mult<-1-vden
  u<-unique(lmm$mrge)
  sel<-which(lmm$mrge!=u[length(u)])
  mrge<-lmm$mrge[sel]
  m2<-length(u)-1
  mult<-1-vden
  # Temperature (air)
  awgts<-aggregate(mult[sel],by=list(mrge),sum)$x
  lwgts<-aggregate(mult[sel],by=list(mrge),length)$x
  wdv <- 0; for (i in 1:m2) wdv<-c(wdv,rep(awgts[i],lwgts[i]))
  wdv<-wdv[-1]
  wgts<-mult[sel]/wdv
  tc2<-aggregate(tc[sel]*wgts,by=list(mrge),sum)$x
  # Temperature (leaf)
  awgts <-aggregate(vden[sel],by=list(mrge),sum)$x
  lwgts <-aggregate(vden[sel],by=list(mrge),length)$x
  wdv<-0; for (i in 1:m2) wdv<-c(wdv,rep(awgts[i],lwgts[i]))
  wdv<-wdv[-1]; wgts<-vden[sel]/wdv
  tleaf2 <- aggregate(tleaf[sel]*wgts,by=list(mrge),sum)$x
  # Conductances (in series)
  gha2<-1/aggregate(1/gha[sel],by=list(mrge),mean)$x  # Need to check
  gt2<-1/aggregate(1/gt,list(lmm$mrge),sum)$x
  # Other variables
  zla2 <- aggregate(zla[sel],list(mrge),mean)$x
  ph2 <-phair(tc2,pk); cp2 <-cpair(tc2)
  mult2<- aggregate(mult[sel],list(mrge),mean)$x
  # PAI needs to be leaf area per metre
  PAI2 <- aggregate(c(PAI,0)/c(lmm$zth,2),list(lmm$mrge),mean)$x
  z2 <- aggregate(z[sel],list(mrge),mean)$x; z2<-c(0,z2,hgt+2)
  # Moisture variables
  Vo2<-aggregate(Vo[sel],list(mrge),mean)$x
  Vflux2<-aggregate(Vflux[sel],list(mrge),mean)$x
  lambda2<- -42.575*tc2+44994
  L2<-aggregate(c(L,0),list(lmm$mrge),mean)$x
  H2<-aggregate(c(H,0),list(lmm$mrge),mean)$x
  TT2 <- aggregate(TT[sel],list(mrge),mean)$x
  return(list(tc=tc2,tleaf=tleaf2,gha=gha2,gt=gt2,zla=zla2,z=z2,ph=ph2,cp=cp2,Vo=Vo2,
              Vflux=Vflux2,lambda=lambda2,L=L2,H=H2,vden=1-mult2,m=m2,PAI=PAI2,TT=TT2))
}
#' Interpolates values from merged canopy layers
#'
#' @description interpolates temperature of vapour concentrations from merged canopy layers
#' to derlive values for the original canopy layers
#'
#' @param y1 Conductance times from ground to each merged canopy node (s)
#' @param y2 Conductance times from ground to each unmerged canopy node (s)
#' @param x1 value to be interpolated (temperature or vapour)
#' @return interpolated value (temperature of vapour)
#' @export
layerinterp <- function(y1, y2, x1) {
  x2 <- 0
  for (i in 1:length(y2)) {
    sel <- which(y1 < y2[i])
    d1 <- abs(y1[max(sel)]-y2[i])
    d2 <- abs(y1[max(sel)+1]-y2[i])
    p1 <- d2/(d1+d2)
    p2 <- d1/(d1+d2)
    x2[i] <- p1*x1[max(sel)]+p2*x1[max(sel)+1]
  }
  x2
}
#' Calculates soil heat conductivity and capacity
#'
#' @description Calculates soil heat conductivity and capacity from soil properites
#' @param timestep model time step (s)
#' @param m number of soil layers
#' @param sdepth depth of deepest soil layer (m)
#' @param theta volumetric soil water fraction (m^3 / m^3)
#' @param frm volumetric soil mineral fraction (m^3 / m^3)
#' @param frq volumetric soil quartz fraction (m^3 / m^3)
#' @param frc mass fraction of clay (kg / kg)
#' @param rho bulk density of soil (Mg / m^3)
#' @return a list with the following components:
#' @return `cd` specific heat capacity of soil (J / m^3 / K)
#' @return `k` thermal conductance of soil (W / m^2 / K)
#' @export
#'
soilk <- function(timestep, m, sdepth = 2, theta = 0.3, frm = 0.3, frq = 0.3, frc = 0.01, rho = 2.65) {
  xx<-(2:(m+1))
  ch<-(2400000*rho/2.64+4180000*theta)
  frs<-frm+frq
  c1<-(0.57+1.73*frq+0.93*frm)/(1-0.74*frq-0.49*frm)-2.8*frs*(1-frs)
  c2<-1.06*rho*theta; c3<-1+2.6*frc^-0.5
  c4<-0.03+0.7*frs^2
  la<-(c1+c2*theta-(c1-c4)*exp(-(c3*theta)^4))
  zmn<-sqrt((la*timestep)/ch)*sqrt(2)
  p<-log(2/zmn,m)
  z<-c(0,0,(sdepth/m^p)*c(1:m)^p)
  cd<-ch*(z[xx+1]-z[xx-1])/(2*timestep)
  k<-la/(z[xx+1]-z[xx])
  return(list(cd=cd, k=k))
}


