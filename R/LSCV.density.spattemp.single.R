LSCV.density.spattemp.single <- function(bands,pp,tt,tlim,sres,tres,grx,gry,grt,kt,inside,xyin,xys,sedge,tedge,parallelise,verbose){
  if(any(bands<=0)) return(NA)
  if(verbose) cat("h =",bands[1],"\b; lambda =",bands[2],"\n")
  h <- bands[1]
  lam <- bands[2]
  
  temp.dens <- kde3d(x=pp$x,y=pp$y,z=tt,h=c(h,h,lam),n=c(sres,sres,tres),lims=c(range(grx),range(gry),kt))
  sq <- matrix(1,sres,sres)
  if(sedge){
    sz <- density.ppp(pp,sigma=h,dimyx=sres,spill=1)
    sq <- sz$edg
    sq[sq>1] <- 1
  }
  sq[!inside] <- NA
  sq <- t(as.matrix(sq))
  tq <- rep(1,tres)
  if(tedge){
    nearedge <- 1:tres
    wellinside <- which(grt>(tlim[1]+4*lam) & grt<(tlim[2]-4*lam))
    if(length(wellinside)>0) nearedge <- nearedge[-wellinside]
    for(i in nearedge) tq[i] <- pnorm(tlim[2],mean=grt[i],sd=lam) - pnorm(tlim[1],mean=grt[i],sd=lam)
  }
  
  if(tedge||sedge){
    for(i in 1:dim(temp.dens$d)[3]) temp.dens$d[,,i] <- temp.dens$d[,,i]/(sq*tq[i])
  }
  # temp.dens <- spattemp.density2(pp,h,lam,tlim,edge,res,outside=NA)
  
  temp.dens.pts <- spattemp.LOO(pp,tt,h,lam,tlim,xyin,xys,sedge,tedge,parallelise=parallelise)
  # temp.dens.pts <- spattemp.density2(pp,h,lam,tlim,edge,res,outside=NA,leaveoneout=TRUE)
  
  return(sum(temp.dens$d^2*xys[1]*xys[2]*(grt[2]-grt[1]),na.rm=TRUE)-2*mean(temp.dens.pts))
}
