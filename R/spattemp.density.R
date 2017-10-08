#' @export
spattemp.density <- function(pp,h=NULL,tt=NULL,lambda=NULL,tlim=NULL,sedge=c("uniform","none"),tedge=sedge,sres=128,tres=NULL,verbose=TRUE){
  if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")
  W <- Window(pp)
  n <- npoints(pp)
  
  sres <- checkit(sres,"'sres'")
  
  sedge <- checkedge(sedge,v=2)
  tedge <- checkedge(tedge,v=2)
  
  WM <- as.mask(W,dimyx=sres)
  inside <- WM$m
  grx <- WM$xcol
  gry <- WM$yrow
  
  if(is.null(tt)) tt <- marks(pp)
  tt <- checktt(tt)
  if(length(tt)!=n) stop(paste("Length of temporal vector does not match number of spatial observations\n   npoints(pp) = ",n,"; length(tt) = ",length(tt),sep=""))
  if(is.null(tlim)) tlim <- range(tt)
  tlim <- checkranin(tlim,tt,"tlim")
  
  if(is.null(tres)){
    tcw <- 1
    kt <- tlim <- c(floor(tlim[1]),ceiling(tlim[2]))
    grt <- tlim[1]:tlim[2]
    tres <- length(grt)
    # kt <- c(tlim[1]-0.5,tlim[2]+0.5)
    # grt <- grt[-which((grt<tlim[1])|(grt>tlim[2]))]
  } else {
    tres <- checkit(tres,"'tres'")
    tcw <- diff(tlim)/tres
    grt <- tlim[1]+0.5*tcw+(0:(tres-1))*tcw
    kt <- c(tlim[1]+0.5*tcw,tlim[2]-0.5*tcw)
  }
  
  if(is.null(h)) h <- OS(pp)
  h <- checkit(h,"'h'")
  if(is.null(lambda)) lambda <- bw.SJ(tt)
  lam <- checkit(lambda,"'lam'")
  
  # fhat <- kde(cbind(pp$x,pp$y,tt),
  #             H=diag(c(h^2,h^2,lam^2)),
  #             xmin=c(min(grx),min(gry),kt[1]),
  #             xmax=c(max(grx),max(gry),kt[2]),
  #             gridsize=c(sres,sres,tres),
  #             supp=4,
  #             verbose=verbose)
  
  if(verbose) message("Calculating trivariate smooth...", appendLF=FALSE)
  fhat <- kde3d(x=pp$x,y=pp$y,z=tt,h=c(h,h,lam),n=c(sres,sres,tres),lims=c(range(grx),range(gry),kt))
  if(verbose) message("Done.")

  if(verbose&&(sedge=="uniform"||tedge=="uniform")) message("Edge-correcting...", appendLF=FALSE)
  sz <- density.ppp(pp,sigma=h,edge=(sedge=="uniform"),dimyx=sres,spill=1)
  sq <- im(matrix(1,sres,sres),xcol=grx,yrow=gry)
  if(sedge=="uniform"){
    sq <- sz$edg
    sq[sq>1] <- 1
  }
  sq[!inside] <- NA
  
  tq <- rep(1,tres)
  if(tedge=="uniform"){
    nearedge <- 1:tres
    wellinside <- which(grt>(tlim[1]+4*lam) & grt<(tlim[2]-4*lam))
    if(length(wellinside)>0) nearedge <- nearedge[-wellinside]
    for(i in nearedge) tq[i] <- pnorm(tlim[2],mean=grt[i],sd=lam) - pnorm(tlim[1],mean=grt[i],sd=lam)
  }
  
  spatial.z <- sz$raw/sq
  spatial.z <- spatial.z/integral(spatial.z)
  temporal.z <- density(tt,bw=lam,from=min(grt),to=max(grt),n=tres)
  temporal.z$y <- temporal.z$y/tq
  if(verbose&&(sedge=="uniform"||tedge=="uniform")) message("Done.")
  
  if(verbose) message("Conditioning on time...", appendLF=FALSE)
  z <- z.cond <- list()
  for(i in 1:tres){
    z[[i]] <- im(t(fhat$d[,,i]),xcol=grx,yrow=gry)
    z[[i]] <- z[[i]]/(sq*tq[i])
    z[[i]][!inside] <- NA
    z.cond[[i]] <- z[[i]]/temporal.z$y[i]
    # z.cond[[i]] <- z.cond[[i]]/integral(z.cond[[i]])
  }
  names(z) <- names(z.cond) <- grt
  if(verbose) message("Done.")
  
  if(sedge=="none") sq <- NULL
  if(tedge=="none") tq <- NULL
  
  final <- list(z=z)
  final$z.cond <- z.cond
  final$h <- h
  final$lambda <- lam
  final$tlim <- tlim #range(grt)
  final$spatial.z <- spatial.z
  final$temporal.z <- temporal.z
  # final$tstep <- tcw
  # final$tbreaks <- seq(tlim[1],tlim[2],length=tres+1)
  # final$tbin <- findInterval(tt,final$tbreaks,all.inside=TRUE)
  final$qs <- sq
  final$qt <- tq
  marks(pp) <- tt
  final$pp <- pp
  final$tgrid <- grt

  class(final) <- "stden"
  return(final)
}
