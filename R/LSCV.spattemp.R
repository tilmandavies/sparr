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



# ## deprecated below (doesn't work, something wrong) ##
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

  