#' Internal function for computing forced convection
.gforced<-function(lw,H,Hmin) {
  d<-0.71*lw
  dT<-0.73673*(d*H^4)^0.2
  dT2<-0.80936*(d*H^4)^0.2
  gha<-0.05*(dT/d)^0.25
  gha2<-0.025*(dT2/d)^0.25
  sel<-which(H<0)
  gha[sel]<-gha2[sel]
  # gmin
  dT<-0.73673*(d*Hmin^4)^0.2
  gmn<-0.05*(dT/d)^0.25
  sel<-which(gha<gmn)
  gha[sel]<-gmn
  gha
}
#' Internal function for computing turbulent convection
.gturb2<-function(u,zu,z1,z0,hgt,PAI) {
  d<-zeroplanedis(hgt,PAI)
  zm<-roughlength(hgt,PAI)
  zh<-0.2*zm
  if (is.na(z0)) z0<-d+zh
  uf<-(0.4*u)/log((zu-d)/zm)
  g<-(0.4*uf*43)/log((z1-d)/(z0-d))
  g
}
#' Internal function for computing below canopy turbulent convection
.gcanopy2<-function(u,zu,z1,z0,hgt,PAI) {
  # Calculate uh
  d<-zeroplanedis(hgt,PAI)
  zm<-roughlength(hgt,PAI)
  uf<-(0.4*u)/log((zu-d)/zm)
  uh<-(uf/0.4)*log((hgt-d)/zm)
  # Calculate g
  a<-attencoef(hgt,PAI,0.2,0.5,1)
  l_m<-mixinglength(hgt,PAI)
  tp<-l_m*0.5*43*uh*a
  e0<-exp(-a*(z0/hgt-1))
  e1<-exp(-a*(z1/hgt-1))
  g<-tp/(e0-e1)
  g
}
#' Internal function to compute canopy and ground radiation absorption
.cgradabs<-function(Rsw,tc,skyem,dp,tme,pai,x,lat,long,alb,galb,slope,
                    aspect,merid,dst,clump) {
  jd<-jday(tme=tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  sa<-solalt(lt,lat,long,jd,merid,dst)
  trb<-cantransdir(pai,x,sa,alb,clump)
  trd<-cantransdif(pai,alb,clump)
  if (is.na(dp[1])) dp<-difprop(Rsw,jd,lt,lat,long,TRUE,TRUE,merid,dst)
  # Ground radiation absorption (sw)
  si<-solarcoef(slope,aspect,lt,lat,long,jd,merid,dst)
  zen<-90-sa
  sim<-si/cos(zen*(pi/180))
  gRsw<-(1-galb)*((1-dp)*trb*sim*Rsw+dp*trd*Rsw)
  # Canopy radiation absorption (sw)
  k<-sqrt((x^2+(tan(zen*(pi/180))^2)))/(x+1.774*(x+1.182)^(-0.733))
  k[sa<=0]<-1
  k[k>10]<-10
  cRb<-(1-dp)*(1-trb)*k*Rsw
  cRd<-dp*trd*Rsw
  cRsw<-(1-alb)*(cRb+cRd)
  # Total
  sb<-5.67*10^-8
  Rem<-0.97*sb*(tc+273.15)^4
  radabs<-gRsw+cRsw+skyem*0.97*Rem
  radabs
}
#' Internal function to compute ground radiation absorption
.gradabs<-function(Rsw,tc,skyem,dp,tme,pai,x,lat,long,alb,galb,slope,
                   aspect,merid,dst,clump) {
  jd<-jday(tme=tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  sa<-solalt(lt,lat,long,jd,merid,dst)
  trb<-cantransdir(pai,x,sa,alb,clump)
  trd<-cantransdif(pai,alb,clump)
  if (is.na(dp[1])) dp<-difprop(Rsw,jd,lt,lat,long,TRUE,TRUE,merid,dst)
  # Ground radiation absorption (sw)
  si<-solarcoef(slope,aspect,lt,lat,long,jd,merid,dst)
  zen<-90-sa
  sim<-si/cos(zen*(pi/180))
  gRsw<-(1-galb)*((1-dp)*trb*sim*Rsw+dp*trd*Rsw)
  gRlw<-canlw(tc,pai,0.03,skyem,clump)$lwabs
  radabs<-gRsw+gRlw
  radabs
}
#' Internal function to compute leaf radiation absorption
.lradabs<-function(Rsw,tc,skyem,dp,tme,paia,x,lat,long,alb,merid,dst,clump) {
  jd<-jday(tme=tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  sa<-solalt(lt,lat,long,jd,merid,dst)
  zen<-90-sa
  trb<-cantransdir(paia,x,sa,alb,clump)
  trd<-cantransdif(paia,alb,clump)
  k<-sqrt((x^2+(tan(zen*(pi/180))^2)))/(x+1.774*(x+1.182)^(-0.733))
  k[sa<=0]<-1
  k[k>10]<-10
  cRb<-(1-dp)*(1-trb)*k*Rsw
  cRd<-dp*trd*Rsw
  cRsw<-(1-alb)*(cRb+cRd)
  cRlw<-canlw(tc,paia,0.03,skyem,clump)$lwabs
  radabs<-cRsw+cRlw
  radabs
}
#' Internal function to calculate leaf stomatal conductance
.gstomatal<-function(Rsw,gsmax) {
  rpar<-Rsw*4.6
  gs<-(gsmax*rpar)/(rpar+100)
  gs
}
#' Internal function to apply Penman-Monteith equation
.PenMont <- function(tc,pk,ea,radabs,gHa,gs,G,Bowen,upper,swet,lims) {
  # Estimate temperature simply using Bowen ratio
  sb<-5.67*10^-8
  Rem<-0.97*sb*(tc+273.15)^4
  Rnet<-radabs-Rem
  cp<-cpair(tc)
  H<-(Bowen*Rnet)/(1+Bowen)
  sel<-which(gs==0)
  tx<-tc+(H/(cp*gHa))
  tx[sel]<-tc[sel]
  es<-satvap(tc)
  tdew<-dewpoint(ea,tc)
  delta<-(4098*es)/(tx+237.3)^2
  gv<-1/(1/gHa+1/gs)
  gHr<-gHa+(4*0.97*sb*(tx+273.15)^3)/cp
  Rem<-0.97*sb*(tc+273.15)^4
  lambda<-(-42.575*tc+44994)
  m<-lambda*swet*(gv/pk)
  T0<-tc+((radabs-Rem-m*(es-ea)-G)/(cp*gHr+m*delta))
  if (lims) {
    sel<-which(T0<tdew)
    T0[sel]<-tdew[sel]
    tmx<-tc+upper
    sel<-which(T0>tmx)
    T0[sel]<-tmx[sel]
  }
  T0
}
#' Internal function to apply Penman-Monteith equation below canopy
.PenBelow<- function(tc,t0,tmx,rh,pk,gtt,gt0,gha,gv,gL,radabs,leafdens,Bowen,lims) {
  cp<-cpair(tc)
  # Air temperature expressed as leaf temperature
  aL<-(gtt*(tc+273.15)+gt0*(t0+273.15))/(gtt+gt0)
  bL<-(leafdens*gL)/(gtt+gt0)
  # Vapour pressures
  es<-satvap(tc)
  eref <- (rh/100)*es
  esoil<-satvap(t0)
  # add a small correction to improve delta estimate
  sb<-5.67*10^-8
  Rnet<-radabs-0.97*sb*(tc+273.15)^4
  H<-(Bowen*Rnet)/(1+Bowen)
  tle<-H/(cp*gha)
  tle[tle< -5]<- -5
  tle[tle> 10]<- 10
  tle<-tle+tc
  wgt1<-ifelse(leafdens<1,leafdens/2,1-0.5/leafdens)
  wgt2<-1-wgt1
  sel<-which(gv==0)
  tx<-wgt1*tle+wgt2*tc
  tx[sel]<-tc[sel]
  delta <- 4098*(satvap(tx))/(tx+237.3)^2
  ae<-(gtt*eref+gt0*esoil+gv*es)/(gtt+gt0+gv)
  be<-(gv*delta)/(gtt+gt0+gv)
  # Sensible heat
  bH<-gha*cp
  # Latent heat
  lambda <- (-42.575*tc+44994)
  aX<-((lambda*gv)/pk)*(es-ae)
  bX<-((lambda*gv)/pk)*(delta-be)
  aX[aX<0]<-0
  bX[bX<0]<-0
  # Emmited radiation
  aR<-sb*0.93*aL^4
  bR<-4*0.97*sb*(aL^3*bL+(tc+273.15)^3)
  # Leaf temperature
  dTL <- (radabs-aR-aX)/(1+bR+bX+bH)
  # tz pass 1
  tn<-aL-273.15+bL*dTL
  tleaf<-tn+dTL
  # new vapour pressure
  eanew<-ae+be*dTL
  eanew[eanew<0.01]<-0.01
  tmn<-dewpoint(eanew,tn,ice = TRUE)
  tmn<-pmax(tmn,tc-7)
  esnew<-satvap(tn)
  sel<-which(eanew>esnew); eanew[sel]<-esnew[sel]
  rh<-(eanew/esnew)*100
  if (lims) {
    # Set both tair and tleaf so as not to drop below dewpoint
    sel<-which(tleaf<tmn); tleaf[sel]<-tmn[sel]
    sel<-which(tn<tmn); tn[sel]<-tmn[sel]
    # Set upper limits
    sel<-which(tleaf>tmx); tleaf[sel]<-tmx[sel]
    sel<-which(tn>tmx); tn[sel]<-tmx[sel]
  }
  return(list(tleaf=tleaf,tn=tn,rh=rh))
}
#' Internal function to calculate the fraction of radiation intercepted by an animal
.Fanim<-function(dlen,dhgt,sa) {
  if (dlen > dhgt) {
    x<-dhgt/dlen
    theta<-sa*(pi/180)
  } else {
    x<-dlen/dhgt
    theta<-(90-sa)*(pi/180)
  }
  if (x == 1) {
    k<-1
  } else {
    top<-sqrt(1+(x^2-1)*cos(theta)*cos(theta))
    asi<-asin(sqrt(1-x^2))
    btm<-2*x+((2*asi)/sqrt(1-x^2))
    k<-top/btm
  }
  k
}
#' Apply height adjustment to wind speed
#'
#' @description Function to apply height adjustment to wind speed,
#' assuming a reference short grass surface
#'
#' @param u wind speed at height `zi` (m/s)
#' @param zi height of input wind speed (m)
#' @param zo height for which wind speed is required (m)
#' @return wind speed at height `zo` (m/s)
#'
#' @examples
#' zo<-c(20:100)/10
#' uo<-windadjust(1,2,zo)
#' plot(zo~uo, type="l", xlab="Wind speed (m/s)", ylab="Height (m)")
#' @export
windadjust <- function(ui,zi,zo) {
  uf<-(0.4*ui)/log((zi-0.08)/0.01476)
  uo<-(uf/0.4)*log((zo-0.08)/0.01476)
  uo
}
#' Compute temperature of heat exchange surface
#'
#' @description Function to compute the temperature of the heat exchange
#' surface of a vegetated canopy.
#'
#' @param climdata a dataframe of hourly weather data formated and with units
#' as per the internal dataset `climdata`
#' @param lat the latitude of the location (decimal degrees)
#' @param long the longitude of the location (decimal degrees)
#' @param hgt canopy height (m)
#' @param pai the total one sided area of canopy elements per unit ground area (see details)
#' @param x leaf distribution angle coefficient
#' @param gsmax maximum stomatal conductance of leaves (mol / m^2 / s)
#' @param slope the slope of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param aspect the aspect of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param alb the albedo of the canopy surface (either the combined ground and
#' canopy albedo if `method = S` or just the canopy albedo)
#' @param galb ground surface albedo
#' @param zu the height above ground of wind speeds in `climdata` (m)
#' @param Bowen Optional parameter specifying the Bowen Ratio of the surface. Used to improve estimates
#' of temperature in application of the Penman-Monteith equation
#' @param upper optional upper limit to temperature offset (difference between reference
#' and canopy surface temperature cannot exceed this value). Ignored if `lims` = FALSE
#' @param umin optional minimum wind speed for computing conductances (avoids conductances being too low)
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT)
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#' @param clump clumpiness factor (0-1, see details)
#' @param method if `S` treats the vegetation surface as a flat surface (see details)
#' @param lims optional logical indicating whether to limit temperatures by `upper` and dewpoint temperature
#' @return a vector of temperatures (deg C) of the canopy heat exchange surface
#'
#' @details if `pai` is unknown it can be estimated as -ln(1-fractional
#' canopy cover). if `clump` = 0 the canopy is assumed entirely uniform
#' and radiation transmission is as for a turbid medium. As `clump`
#' approaches 1, the canopy is assumed to be increasingly patchy, such
#' that a greater proportion of reaches the ground without being obscured
#' by leaves. If `method` = S, the radiation intercepted by the canopy is assumed to
#' be that for a flat surface. If `method` is not S, the radiation absorption
#' by the gorund and canopy surface are computed separately accounting for the
#' inclination of the ground surface and the distribution of leaf angles.
#'
#' @export
#'
#' @examples
#' # Compute temperature a surface
#' Tref<-Thes(climdata,50.2178,-5.32656,hgt=0.5,pai=2,x=1,gsmax=0.33,alb=0.23,method="C",lims=TRUE)
#' plot(Tref,type="l")
Thes<-function(climdata,lat,long,hgt,pai,x=1,gsmax=0.33,slope=0,aspect=0,alb=0.23,galb=0.15,
               zu=2,Bowen=1.5,upper=25,umin=0.5,merid=0,dst=0,clump=0,method="S",lims="FALSE") {
  if(zu<hgt) stop("cannot compute when zu<hgt - apply windadjust to derive wind above canopy")
  # extract data
  tc<-climdata$temp
  rh<-climdata$relhum
  u<-climdata$windspeed
  Rsw<-climdata$swrad
  skyem<-climdata$skyem
  pk<-climdata$pres
  tme<-as.POSIXlt(climdata$obs_time,tz="GMT")
  dp<-climdata$difrad/climdata$swrad
  dp[is.na(dp)]<-0.5
  dp[dp<0]<-0
  dp[dp>1]<-0
  # Calculate
  ea<-satvap(tc)*(rh/100)
  u[u<umin]<-umin
  gHa<-.gturb2(u,zu,zu,NA,hgt,pai)
  gs<-3*.gstomatal(Rsw,gsmax)
  sb<-5.67*10^-8
  Rem<-0.97*sb*(tc+273.15)^4
  Rnet<-Rsw+(1-skyem)*Rem
  sel<-which(Rnet>0)
  G<-0.5*Rnet
  G[sel]<-0.1*Rnet[sel]
  if (method=="S") {
    radabs<-(1-alb)*Rsw+skyem*0.97*Rem
  } else radabs<-.cgradabs(Rsw,tc,skyem,dp,tme,pai,x,lat,long,alb,galb,slope,aspect,merid,dst,clump)
  Th<-.PenMont(tc,pk,ea,radabs,gHa,gs,G,Bowen,upper,1,lims)
  Th
}
#' Compute temperature above canopy
#'
#' @description Function to compute air temperature above canopy.
#'
#' @param climdata a dataframe of hourly weather data formated and with units
#' as per the internal dataset `climdata`
#' @param z height above canopy for which temperature estimate is required (m)
#' @param lat the latitude of the location (decimal degrees)
#' @param long the longitude of the location (decimal degrees)
#' @param hgt canopy height (m)
#' @param pai the total one sided area of canopy elements per unit ground area (see details)
#' @param x leaf distribution angle coefficient
#' @param gsmax maximum stomatal conductance of leaves (mol / m^2 / s)
#' @param slope the slope of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param aspect the aspect of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param alb the albedo of the canopy surface (either the combined ground and
#' canopy albedo if `method = S` or just the canopy albedo)
#' @param galb ground surface albedo
#' @param zu the height above ground of wind speeds in `climdata` (m)
#' @param Bowen Optional parameter specifying the Bowen Ratio of the surface. Used to improve estimates
#' of temperature in application of the Penman-Monteith equation
#' @param upper optional upper limit to temperature offset (difference between reference
#' and canopy surface temperature cannot exceed this value). Ignored if `lims` = FALSE
#' @param umin optional minimum wind speed for computing conductances (avoids conductances being too low)
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT)
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#' @param clump clumpiness factor (0-1, see details)
#' @param method if `S` treats the vegetation surface as a flat surface (see details)
#' @param lims optional logical indicating whether to limit temperatures by `upper` and dewpoint temperature
#' @return a vector of air temperatures (deg C) above canopy
#'
#' @details if `pai` is unknown it can be estimated as -ln(1-fractional
#' canopy cover). if `clump` = 0 the canopy is assumed entirely uniform
#' and radiation transmission is as for a turbid medium. As `clump`
#' approaches 1, the canopy is assumed to be increasingly patchy, such
#' that a greater proportion of reaches the ground without being obscured
#' by leaves. If `method` = S, the radiation intercepted by the canopy is assumed to
#' be that for a flat surface. If `method` is not S, the radiation absorption
#' by the gorund and canopy surface are computed separately accounting for the
#' inclination of the ground surface and the distribution of leaf angles.
#'
#' @export
#'
#' @examples
#' # Compute temperature offset of reference surface
#' Tref<-Tabove(climdata,1,50.2178,-5.32656,hgt=0.5,pai=2,gsmax=0.33,alb=0.23)
#' dTRef<-Tref-climdata$temp
#' # Compute temperature offset of different surface
#' Th<-Tabove(climdata,1,50.2178,-5.32656,hgt=0.25,pai=1.5,gsmax=0.43,alb=0.15)
#' dTh<-Th-climdata$temp
#' # Check whether offsets are linearly related to one another
#' plot(dTh~dTref,pch=15)
#' abline(lm(dTh~dTref),lwd=2,col="red")
Tabove<-function(climdata,z,lat,long,hgt,pai,x=1,gsmax=0.33,slope=0,aspect=0,alb=0.23,galb=0.15,
                 zu=2,Bowen=1.5,upper=25,umin=0.5,merid=0,dst=0,clump=0,method="S",lims="FALSE") {
  # extract data
  tc<-climdata$temp
  pk<-climdata$pres
  u<-climdata$windspeed
  u[u<umin]<-umin
  gHa<-.gturb2(u,zu,zu,NA,hgt,pai)
  ph<-phair(tc,pk)
  d<-zeroplanedis(hgt,pai)
  zm<-roughlength(hgt,pai)
  Th<-Thes(climdata,lat,long,hgt,pai,x,gsmax,slope,aspect,alb,galb,zu,Bowen,upper,umin,merid,dst,clump,method,lims)
  Tz<-Th-(gHa*(Th-tc)/(0.4*ph))*log((z-d)/zm)
  Tz
}
#' Compute ground surface temperature
#'
#' @description Function to compute the ground surface temperature using the Penman-Monteith equation
#'
#' @param climdata a dataframe of hourly weather data formated and with units
#' as per the internal dataset `climdata`
#' @param lat the latitude of the location (decimal degrees)
#' @param long the longitude of the location (decimal degrees)
#' @param hgt canopy height (m)
#' @param pai the total one sided area of canopy elements per unit ground area (see details)
#' @param x leaf distribution angle coefficient
#' @param gsmax maximum stomatal conductance of leaves (mol / m^2 / s)
#' @param slope the slope of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param aspect the aspect of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param alb the albedo of the canopy surface (either the combined ground and
#' canopy albedo if `method = S` or just the canopy albedo)
#' @param galb ground surface albedo
#' @param zu the height above ground of wind speeds in `climdata` (m)
#' @param Bowen Optional parameter specifying the Bowen Ratio of the surface. Used to improve estimates
#' of temperature in application of the Penman-Monteith equation
#' @param upper optional upper limit to temperature offset (difference between reference
#' and canopy surface temperature cannot exceed this value). Ignored if `lims` = FALSE
#' @param umin optional minimum wind speed for computing conductances (avoids conductances being too low)
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT)
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#' @param clump clumpiness factor (0-1, see details)
#' @param method if `S` treats the vegetation surface as a flat surface (see details)
#' @param lims optional logical indicating whether to limit temperatures by `upper` and dewpoint temperature
#' @return a vector of ground surface temperatures (deg C)
#'
#' @details if `pai` is unknown it can be estimated as -ln(1-fractional
#' canopy cover). if `clump` = 0 the canopy is assumed entirely uniform
#' and radiation transmission is as for a turbid medium. As `clump`
#' approaches 1, the canopy is assumed to be increasingly patchy, such
#' that a greater proportion of reaches the ground without being obscured
#' by leaves. If `method` = S, the radiation intercepted by the canopy is assumed to
#' be that for a flat surface. If `method` is not S, the radiation absorption
#' by the ground and canopy surface are computed separately accounting for the
#' inclination of the ground surface and the distribution of leaf angles.
#'
#' @export
#'
#' @examples
#' # Compute temperature a surface
#' Tg<-Tground(climdata,50.2178,-5.32656,hgt=0.5,pai=2,x=1,gsmax=0.33,alb=0.23)
#' plot(Tg,type="l")
Tground<-function(climdata,lat,long,hgt,pai,x=1,gsmax=0.33,slope=0,aspect=0,alb=0.23,galb=0.15,
                  zu=2,Bowen=1.5,upper=25,umin=0.5,merid=0,dst=0,clump=0,method="S",lims=TRUE) {
  if(zu<hgt) stop("cannot compute when zu<hgt - apply windadjust to derive wind above canopy")
  # extract data
  tc<-climdata$temp
  rh<-climdata$relhum
  u<-climdata$windspeed
  Rsw<-climdata$swrad
  skyem<-climdata$skyem
  pk<-climdata$pres
  tme<-as.POSIXlt(climdata$obs_time,tz="GMT")
  dp<-climdata$difrad/climdata$swrad
  dp[is.na(dp)]<-0.5
  dp[dp<0]<-0
  dp[dp>1]<-0
  # Calculate non-conductivities
  ea<-satvap(tc)*(rh/100)
  u[u<umin]<-umin
  d<-zeroplanedis(hgt,pai)
  zm<-roughlength(hgt,pai)
  uf<-(0.4*u)/log((zu-d)/zm)
  uh<-(uf/0.4)*log((hgt-d)/zm)
  # Calculate conductivities
  gh2<-.gturb2(u,zu,zu,hgt,hgt,pai)
  g0h<-.gcanopy2(u,zu,hgt,0,hgt,pai)
  gHa<-1/(1/gh2+1/g0h)
  # Calculate G
  sb<-5.67*10^-8
  Rem<-0.97*sb*(tc+273.15)^4
  Rnet<-Rsw+(1-skyem)*Rem
  sel<-which(Rnet>0)
  G<-0.25*Rnet
  G[sel]<-0.1*Rnet[sel]
  # Calculate radiation asborbption
  if (method=="S") {
    radabs<-(1-alb)*Rsw*exp(-pai)+skyem*0.97*Rem
  } else radabs<-.gradabs(Rsw,tc,skyem,dp,tme,pai,x,lat,long,alb,galb,slope,aspect,merid,dst,clump)
  Tg<-.PenMont(tc,pk,ea,radabs,gHa,gHa,G,Bowen,upper,1,lims)
  Tg
}
#' Compute below-canopy air temperature
#'
#' @description Function to compute below canopy air temperatures
#'
#' @param climdata a dataframe of hourly weather data formated and with units
#' as per the internal dataset `climdata`
#' @param z the height (below canopy) for which temperature estimates are required (m)
#' @param lat the latitude of the location (decimal degrees)
#' @param long the longitude of the location (decimal degrees)
#' @param hgt canopy height (m)
#' @param pai the total one sided area of canopy elements per unit ground area (see details)
#' @param paia optionally, the total one sided area of canopy elements per unit ground area above z (see details)
#' @param x leaf distribution angle coefficient
#' @param gsmax maximum stomatal conductance of leaves (mol / m^2 / s)
#' @param lw average leaf width (m)
#' @param slope the slope of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param aspect the aspect of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param alb the albedo of the canopy surface (either the combined ground and
#' canopy albedo if `method = S` or just the canopy albedo)
#' @param galb ground surface albedo
#' @param zu the height above ground of wind speeds in `climdata` (m)
#' @param Bowen Optional parameter specifying the Bowen Ratio of the surface. Used to improve estimates
#' of temperature in application of the Penman-Monteith equation
#' @param upper optional upper limit to temperature offset (difference between reference
#' and canopy surface temperature cannot exceed this value). Ignored if `lims` = FALSE
#' @param umin optional minimum wind speed for computing conductances (avoids conductances being too low)
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT)
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#' @param clump clumpiness factor (0-1, see details)
#' @param method if `S` treats the vegetation surface as a flat surface (see details)
#' @param lims optional logical indicating whether to limit temperatures by `upper` and dewpoint temperature
#' @return a list of the following elements:
#' @return `tleaf` a vector of leaf surface temperatures (deg C)
#' @return `tn` a vector of air temperatures (deg C)
#' @return `rh` a vector of relative humidities (percentage)
#'
#' @details if `pai` is unknown it can be estimated as -ln(1-fractional
#' canopy cover). If `paia` is unspecified vertically uniform foliage distribution is assumed.
#' if `clump` = 0 the canopy is assumed entirely uniform
#' and radiation transmission is as for a turbid medium. As `clump`
#' approaches 1, the canopy is assumed to be increasingly patchy, such
#' that a greater proportion of reaches the ground without being obscured
#' by leaves. If `method` = S, the radiation intercepted by the canopy is assumed to
#' be that for a flat surface. If `method` is not S, the radiation absorption
#' by the ground and canopy surface are computed separately accounting for the
#' inclination of the ground surface and the distribution of leaf angles.
#'
#' @export
#'
#' @examples
#' # Compute temperatures below canopy
#' Tln<-Tcanopy(climdata,0.25,50.2178,-5.32656,hgt=0.5,pai=2,x=1,gsmax=0.33,alb=0.23)
#' plot(Tln$tn,type="l") # air temperature
#' plot(Tln$tleaf,type="l") # leaf temperature
Tcanopy<-function(climdata,z,lat,long,hgt,pai,paia=NA,x=1,gsmax=0.33,lw=0.05,slope=0,aspect=0,alb=0.23,galb=0.15,
                  zu=2,Bowen=1.5,upper=25,umin=0.5,merid=0,dst=0,clump=0,method="C",lims=TRUE) {
  if(zu<hgt) stop("cannot compute when zu<hgt - apply windadjust to derive wind above canopy")
  # extract data
  tc<-climdata$temp
  rh<-climdata$relhum
  u<-climdata$windspeed
  Rsw<-climdata$swrad
  skyem<-climdata$skyem
  pk<-climdata$pres
  tme<-as.POSIXlt(climdata$obs_time,tz="GMT")
  dp<-climdata$difrad/climdata$swrad
  dp[is.na(dp)]<-0.5
  dp[dp<0]<-0
  dp[dp>1]<-0
  # Calculate non-conductivities
  ea<-satvap(tc)*(rh/100)
  u[u<umin]<-umin
  d<-zeroplanedis(hgt,pai)
  zm<-roughlength(hgt,pai,zm0=hgt/100)
  uf<-(0.4*u)/log((zu-d)/zm)
  uh<-(uf/0.4)*log((hgt-d)/zm)
  # Calculate wind speed and radiation below canopy
  a<-attencoef(hgt,pai)
  uz<-uh*exp(a*((z/hgt)-1))
  dp<-climdata$difrad/climdata$swrad
  dp[is.na(dp)]<-0.5
  dp[dp<0]<-0
  dp[dp>1]<-0
  jd<-jday(tme=tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  # compute paia
  if (is.na(paia[1])) paia<-((hgt-z)/hgt)*pai
  leafdens<-paia/(hgt-z)
  rsw<-cansw(Rsw,dp,jd,lt,lat,long,paia,x,alb,TRUE,TRUE,merid,dst,clump=clump)
  th<-Tabove(climdata,hgt,lat,long,hgt,pai,x,gsmax,slope,aspect,alb,galb,zu,Bowen,upper,umin,
             merid,dst,clump,method,lims)
  tmx<-Thes(climdata,lat,long,hgt,pai,x,gsmax,slope,aspect,alb,galb,zu,Bowen,upper,umin,
            merid,dst,clump,method,lims)
  t0<-Tground(climdata,lat,long,hgt,pai,x,gsmax,slope,aspect,alb,galb,
              zu,Bowen,upper,umin,merid,dst,clump,method,lims)
  # Set minimum ground tmep
  ea<-satvap(tc)*(rh/100)
  tmn<-dewpoint(ea,tc)
  sel<-which(t0<tmn)
  t0[sel]<-tmn[sel]
  radabs<-.lradabs(Rsw,tc,skyem,dp,tme,paia,x,lat,long,alb,merid,dst,clump)
  # Calculate conductivities
  gtt<-.gcanopy2(u,zu,hgt,z,hgt,pai)
  gt0<-.gcanopy2(u,zu,z,0,hgt,pai)
  gha<-1.4*0.135*sqrt(uz/(0.71*lw))
  # Calculate min conductivity
  sb<-5.67*10^-8
  rnet<-radabs-sb*0.97*(tc+273.15)^4
  H<-(Bowen*rnet)/(1+Bowen)
  gmn<-.gforced(lw,H,5)
  sel<-which(gha<gmn)
  gha[sel]<-gmn[sel]
  gs<-.gstomatal(rsw,gsmax)
  gv<-1/(1/gha+1/gs)
  ph<-phair(tc,pk)
  gL<-1/(1/gha+1/(uz*ph))
  # Calculate temperatures
  Tzl<-.PenBelow(tc,t0,tmx,rh,pk,gtt,gt0,gha,gv,gL,radabs,leafdens,Bowen,lims)
  Tzl
}
#' Compute ectotherm body temperature
#'
#' @description Function to compute ectotherm body temperature
#'
#' @param climdata a dataframe of hourly weather data formated and with units
#' as per the internal dataset `climdata`
#' @param z the height (below canopy) for which temperature estimates are required (m)
#' @param dlen length of animal (m)
#' @param dhgt height of animal (m)
#' @param lat the latitude of the location (decimal degrees)
#' @param long the longitude of the location (decimal degrees)
#' @param hgt canopy height (m)
#' @param pai the total one sided area of canopy elements per unit ground area (see details)
#' @param paia optionally, the total one sided area of canopy elements per unit ground area above z (see details)
#' @param x leaf distribution angle coefficient
#' @param gsmax maximum stomatal conductance of leaves (mol / m^2 / s)
#' @param lw average leaf width (m)
#' @param slope the slope of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param aspect the aspect of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param alb the albedo of the canopy surface (either the combined ground and
#' canopy albedo if `method = S` or just the canopy albedo)
#' @param galb ground surface albedo
#' @param dalb albedo of animal
#' @param skinwet proportion of animal surface acting like a wet surface
#' @param zu the height above ground of wind speeds in `climdata` (m)
#' @param Bowen Optional parameter specifying the Bowen Ratio of the surface. Used to improve estimates
#' of temperature in application of the Penman-Monteith equation
#' @param upper optional upper limit to temperature offset (difference between reference
#' and canopy surface temperature cannot exceed this value). Ignored if `lims` = FALSE
#' @param umin optional minimum wind speed for computing conductances (avoids conductances being too low)
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT)
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#' @param clump clumpiness factor (0-1, see details)
#' @param method if `S` treats the vegetation surface as a flat surface (see details)
#' @param lims optional logical indicating whether to limit temperatures by `upper` and dewpoint temperature
#' @return a vector of animal body temperatures
#'
#' @details if `pai` is unknown it can be estimated as -ln(1-fractional
#' canopy cover). If `paia` is unspecified vertically uniform foliage distribution is assumed.
#' if `clump` = 0 the canopy is assumed entirely uniform
#' and radiation transmission is as for a turbid medium. As `clump`
#' approaches 1, the canopy is assumed to be increasingly patchy, such
#' that a greater proportion of reaches the ground without being obscured
#' by leaves. If `method` = S, the radiation intercepted by the canopy is assumed to
#' be that for a flat surface. If `method` is not S, the radiation absorption
#' by the ground and canopy surface are computed separately accounting for the
#' inclination of the ground surface and the distribution of leaf angles.
#'
#' @export
#'
#' @examples
#' # Compute body temperatures of a toad near the surface below canopy
#' Ttoad<-Tbody(climdata,0.05,0.08,0.06,50.2178,-5.32656,0.5,2)
#' plot(Ttoad,type="l")
Tbody<-function(climdata,z,dlen,dhgt,lat,long,hgt,pai,paia=NA,x=1,gsmax=0.33,lw=0.05,
                slope=0,aspect=0,alb=0.23,galb=0.15,dalb=0.23,skinwet=1,
                zu=2,Bowen=1.5,upper=25,umin=0.5,merid=0,dst=0,clump=0,method="C",lims=TRUE) {
  trd<-1
  trb<-1
  Rlw<-canlw(climdata$temp,0,0.03,climdata$skyem,clump)$lwabs
  ea<-(climdata$relhum/100)*satvap(climdata$temp)
  tme<-as.POSIXlt(climdata$obs_time,tz="UTC")
  jd<-jday(tme=tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  sa<-solalt(lt,lat,long,jd,merid,dst)
  # Work out air temperatures and canopy transmission
  if (z>hgt) {
    ta<-Tabove(climdata,z,lat,long,hgt,pai,x,gsmax,slope,aspect,alb,galb,zu,Bowen,
               upper,umin,merid,dst,clump,method,lims)

  } else if (z==hgt) {
    ta<-Thes(climdata,lat,long,hgt,pai,x,gsmax,slope,aspect,alb,galb,zu,Bowen,
             upper,umin,merid,dst,clump,method,lims)
  } else {
    tzn<-Tcanopy(climdata,z,lat,long,hgt,pai,paia,x,gsmax,lw,slope,aspect,alb,galb,
                 zu,Bowen,upper,umin,merid,dst,clump,method,lims)
    ta<-tzn$tn
    rh<-tzn$rh
    ea<-(rh/100)*satvap(ta)
    trb<-cantransdir(pai,x,sa,alb,clump)
    trd<-cantransdif(pai,alb,clump)
    Rlw<-canlw(climdata$temp,pai,0.03,climdata$skyem,clump)$lwabs
  }
  # work out fraction of beam radiation absorbed (assuming a spheroid)
  Fa<-.Fanim(dlen,dhgt,sa)
  # Calculate radiation absorbed
  Rd<-climdata$difrad
  Rb<-climdata$swrad-Rd
  ze<-(90-sa)*(pi/180)
  Rb<-(climdata$swrad-Rd)/cos(ze)
  Rb[Rb<0]<-0
  Rb[Rb>1352]<-1352
  Rabs<-(1-dalb)*(trd*Rd+trb*Rb*Fa)+Rlw
  # Calculate wind speed at height z
  u<-climdata$windspeed
  u[u<umin]<-umin
  d<-zeroplanedis(hgt,pai)
  zm<-roughlength(hgt,pai,zm0=hgt/100)
  uf<-(0.4*u)/log((zu-d)/zm)
  uh<-(uf/0.4)*log((hgt-d)/zm)
  a<-attencoef(hgt,pai)
  uz<-uh*exp(a*((z/hgt)-1))
  # Calculate boundary layer conductance
  if (dlen > dhgt) {
    volume<-(4/3)*pi*dhgt^2*dlen
  } else volume<-(4/3)*pi*dlen^2*dhgt
  d<-volume^(1/3)
  sb<-5.67*10^-8
  Rnet<-Rabs-sb*0.97*(climdata$temp+273.15)^4
  H<-(Bowen*Rnet)/(1+Bowen)
  gmn<-.gforced(d/0.71,H,5)
  gha<-240*uz^0.6*d^(-0.4)
  sel<-which(gha<gmn)
  gha[sel]<-gmn[sel]
  # Calculate body temp
  pk<-climdata$pres
  Tb<-.PenMont(ta,pk,ea,Rabs,gha,1e10,0,Bowen,upper,swet=skinwet,lims)
  Tb
}
#' Compute terrain influencing factor
#'
#' @description Function to compute terrain influencing factor
#'
#' @param climdata a dataframe of hourly weather data formated and with units
#' as per the internal dataset `climdata`
#' @param lat the latitude of the location (decimal degrees)
#' @param long the longitude of the location (decimal degrees)
#' @param pai the total one sided area of canopy elements per unit ground area (see details)
#' @param slope the slope of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param aspect the aspect of the underlying ground surface (decimal degrees). Ignored if `method` = S.
#' @param merid an optional numeric value representing the longitude (decimal degrees) of the local time zone meridian (0 for GMT)
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#' @return a single numeric value giving a terrain influencing factor
#'
#' @export
#'
#' @examples
#' # Compute influencing factor of an inclined surface
#' iH<-terraincoef(climatedata,50.2178,-5.32656,1,30,180)
#' # Compute influencing factor of a flat surface
#' iRef<-terraincoef(climatedata,50.2178,-5.32656,1,0,0)
#' # Compute ratio
#' iR<-iH/iRef
#' iR
terraincoef<-function(climatedata,lat,long,pai,slope,aspect,merid=0,dst=0) {
  # ~~ compute Smax
  tme<-as.POSIXlt(climdata$obs_time,tz="UTC")
  jd<-jday(tme=tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  sa<-solalt(lt,lat,long,jd,merid,dst)
  zen<-90-sa
  sel<-which(zen==min(zen))
  Smax<-solarcoef(slope,aspect,lt[sel],lat,long,jd[sel],merid,dst)
  # ~~ compute canopy transmission
  tr<-exp(-pai)
  # ~~ compute daytime radiation ratio
  # compute diffuse proportion
  dp<-climdata$difrad/climdata$swrad
  dp[is.na(dp)]<-0.5
  dp[dp>1]<-1
  dp[dp<0]<-0
  # compute daytime ratio
  seld<-which(climdata$swrad>0)
  mdp<-mean(dp[seld])
  rat<-(1-mdp)
  # ~~ compute terrain influencing factor
  ifa<-rat*tr*Smax
  ifa
}
