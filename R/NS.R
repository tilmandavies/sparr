NS <- function(pp, nstar = c("npoints", "geometric"), scaler = c("silverman", "IQR", "sd", "var")){
  if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")
  
  nstar <- processnstar(nstar,pp)
  scaler <- processscaler(scaler,pp)
  
  return(scaler*nstar^(-1/6)) #'was: 2*nstar in denominator prior to adjust
}
  