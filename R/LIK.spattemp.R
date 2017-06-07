#' @export
LIK.spattemp <- function(pp,tt=NULL,tlim=NULL,sedge=c("uniform","none"),tedge=sedge,parallelise=NA,start=NULL,verbose=TRUE){
  if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")
  W <- Window(pp)
  n <- npoints(pp)
  sedge <- checkedge(sedge,v=2)
  tedge <- checkedge(tedge,v=2)
  WM <- as.mask(W,dimyx=64)
  
  if(is.null(tt)) tt <- marks(pp)
  tt <- checktt(tt)
  if(length(tt)!=n) stop(paste("Length of temporal vector does not match number of spatial observations\n   npoints(pp) = ",n,"; length(tt) = ",length(tt),sep=""))
  if(is.null(tlim)) tlim <- range(tt)
  tlim <- checkranin(tlim,tt,"tlim")
  
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

  result <- optim(start,LIK.density.spattemp.single,
                  pp=pp,tt=tt,tlim=tlim,xyin=evalxy.in,
                  xys=c(WM$xstep,WM$ystep),
                  sedge=(sedge=="uniform"),
                  tedge=(tedge=="uniform"),
                  parallelise=parallelise,
                  verbose=verbose)$par
  
  names(result) <- c("h","lambda")
  return(result)
}
