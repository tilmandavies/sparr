#' @rdname NS
#' @export
NS.spattemp <- function(pp, tt = NULL, nstar = "npoints", scaler = c("silverman", "IQR", "sd", "var")){
  if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")
  
  if(is.null(tt)) tt <- marks(pp)
  tt <- checktt(tt)
  if(length(tt)!=npoints(pp)) stop(paste("Length of temporal vector does not match number of spatial observations\n   npoints(pp) = ",npoints(pp),"; length(tt) = ",length(tt),sep=""))
  
  scaler <- processscaler.st(scaler,pp,tt)
  nstar <- processnstar.st(nstar,pp)
  
  result <- scaler*c(nstar^(-1/6),0.9*nstar^(-1/5))
  
  names(result) <- c("h","lambda")
  return(result)
}
