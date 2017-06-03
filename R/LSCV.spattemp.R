#' @export
LSCV.spattemp <- function(pp,tt=NULL,tlim=NULL,sedge=c("uniform","none"),tedge=sedge,sres=64,tres=sres,parallelise=NA,start=NULL,verbose=TRUE){
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
  } else {
    tres <- checkit(tres,"'tres'")
    tcw <- diff(tlim)/tres
    grt <- tlim[1]+0.5*tcw+(0:(tres-1))*tcw
    kt <- c(tlim[1]+0.5*tcw,tlim[2]-0.5*tcw)
  }
  
  if(!is.na(parallelise)){
    if(!is.numeric(parallelise)) stop("'parallelise' must be numeric")
    if(is.null(parallelise)) parallelise <- NA
    parallelise <- round(parallelise[1])
  }
  
  evalxy <- as.matrix(expand.grid(WM$xcol,WM$yrow))
  notin <- !inside.owin(x=evalxy[,1],y=evalxy[,2],w=W)
  evalxy.in <- evalxy[!notin,]
    
  if(is.null(start)) start <- c(OS(pp),bw.SJ(tt))
  if(any(start<0)) stop("invalid starting values in 'start'")
  
  result <- optim(start,LSCV.density.spattemp.single,
                  pp=pp,tt=tt,tlim=tlim,sres=sres,
                  tres=tres,grx=grx,gry=gry,grt=grt,
                  kt=kt,inside=inside,xyin=evalxy.in,
                  xys=c(WM$xstep,WM$ystep),
                  sedge=(sedge=="uniform"),
                  tedge=(tedge=="uniform"),
                  parallelise=parallelise,
                  verbose=verbose)$par

  names(result) <- c("h","lambda")
  return(result)
}


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


spattemp.LOO <- function(pp,tt,h,lambda,tlim,xyin,xys,sedge,tedge,parallelise){
  W <- Window(pp)
  n <- npoints(pp)
  ppmat <- cbind(pp$x,pp$y)
  qs <- qt <- 1
  
  if(is.na(parallelise)){
    loo <- rep(NA,n)
    for(i in 1:n){
      ppt.i <- pp[i]
      ppt.mi <- pp[-i]
      t.i <- tt[i]
      t.mi <- tt[-i]
      
      if(sedge){
        pxy <- kernel2d(xyin[,1]-ppt.i$x, xyin[,2]-ppt.i$y, h)
        qs <- dintegral(pxy,xys[1],xys[2])
      }
      if(tedge) qt <- pnorm(tlim[2],t.i,lambda) - pnorm(tlim[1],t.i,lambda)
      
      ut <- (t.i-t.mi)/lambda
      ivals <- kernel2d(ppt.i$x-ppt.mi$x,ppt.i$y-ppt.mi$y,h)*lambda^(-1)*exp(-0.5*ut^2)/sqrt(2*pi)
      loo[i] <- mean(ivals)/(qs*qt)
    }
  } else {
    registerDoParallel(cores=parallelise)
    loo <- foreach(i=1:n,.packages="spatstat",.combine=c) %dopar% {
      ppt.i <- pp[i]
      ppt.mi <- pp[-i]
      t.i <- tt[i]
      t.mi <- tt[-i]
      
      if(sedge){
        pxy <- kernel2d(xyin[,1]-ppt.i$x, xyin[,2]-ppt.i$y, h)
        qs <- dintegral(pxy,xys[1],xys[2])
      }
      if(tedge) qt <- pnorm(tlim[2],t.i,lambda) - pnorm(tlim[1],t.i,lambda)
      
      ut <- (t.i-t.mi)/lambda
      ivals <- kernel2d(ppt.i$x-ppt.mi$x,ppt.i$y-ppt.mi$y,h)*lambda^(-1)*exp(-0.5*ut^2)/sqrt(2*pi)
      return(mean(ivals)/(qs*qt))
    }
  }
  return(loo)
}




# ## doesn't work below something wrong ##
# LSCV.spattemp.single <- function(hv,pp,tt,tlim,sedge,tedge,sres,tres,para,verbose){
#   if(any(hv<0)) return(NA)
#   if(verbose) cat("h =",hv[1],"\b; lambda =",hv[2],"\n")
#   sttemp <- spattemp.density(pp=pp,h=hv[1],tt=tt,lambda=hv[2],tlim=tlim,sedge=sedge,tedge=tedge,sres=sres,tres=tres,verbose=FALSE)
#   lala <<- sttemp
#   # stslices <- spattemp.slice(sttemp,tt=tt,checkargs=FALSE)
#   # print(range(sttemp$tt))
#   
#   qs <- qt <- rep(1,pp$n)
#   if(tedge=="uniform") qt <- approx(sttemp$tt,sttemp$qt,xout=tt,rule=2)$y
#   if(sedge=="uniform") qs <- safelookup(sttemp$qs,pp,warn=FALSE)
#   
#   # print("lala")
#   
#   sti <- rep(NA,pp$n)
#   if(is.na(para)){
#     #for(i in 1:pp$n) sti[i] <- stsurfs$z[[i]][si$row[i],si$col[i]] - (dnorm(0,sd=hv[1])^2*dnorm(0,sd=hv[2]))/(sttemp$qs[si$row[i],si$col[i]]
#     for(i in 1:pp$n){
#       # print(spattemp.slice(sttemp,tt=tt[i],checkargs=T)$z)
#       # print(safelookup(spattemp.slice(sttemp,tt=tt[i],checkargs=FALSE)$z[[1]],pp[i],warn=FALSE))
#       # print((dnorm(0,sd=hv[1])^2*dnorm(0,sd=hv[2]))/(qs[i]*qt[i]))
#       # 
#       # print(tt[i])
#       sti[i] <- safelookup(spattemp.slice(sttemp,tt=tt[i],checkargs=FALSE)$z[[1]],pp[i],warn=FALSE) - (dnorm(0,sd=hv[1])^2*dnorm(0,sd=hv[2]))/(qs[i]*qt[i])
#     }  
#   } else {
#     registerDoParallel(cores=para)
#     sti <- foreach(i=1:pp$n,.packages="spatstat",.combine=c) %dopar% {
#       return(safelookup(spattemp.slice(sttemp,tt=tt[i],checkargs=FALSE)$z[[1]],pp[i],warn=FALSE) - (dnorm(0,sd=hv[1])^2*dnorm(0,sd=hv[2]))/(qs[i]*qt[i]))
#     }
#   }
#   
#   stsq <- lapply(sttemp$z,function(x) x^2)
#   stint <- sum(Reduce("+",stsq))*sttemp$z[[1]]$xstep*sttemp$z[[1]]$ystep*(sttemp$tt[2]-sttemp$tt[1])
#   
#   result <- stint-2*mean(sti)
#   return(result)
# }

  