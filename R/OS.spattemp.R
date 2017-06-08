#' @rdname OS
#' @export
OS.spattemp <- function(pp, tt = NULL, nstar = "npoints", scaler = c("silverman", "IQR", "sd", "var")){
  if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")
  
  if(is.null(tt)) tt <- marks(pp)
  tt <- checktt(tt)
  if(length(tt)!=npoints(pp)) stop(paste("Length of temporal vector does not match number of spatial observations\n   npoints(pp) = ",npoints(pp),"; length(tt) = ",length(tt),sep=""))
  
  scaler <- processscaler.st(scaler,pp,tt)
  nstar <- processnstar.st(nstar,pp)
  
  RK <- c(1/(4*pi),1/(2*sqrt(pi)))
  d <- 2:1
  V <- (16*gamma((d+8)/2)*d*(d+2))/((d+8)^((d+6)/2)*pi^(d/2))
  result <- scaler*(((RK*d)/(nstar*V))^(1/(d+4)))
  
  names(result) <- c("h","lambda")
  return(result)
}