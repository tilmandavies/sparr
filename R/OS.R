OS <- function(pp, nstar = c("npoints", "geometric"), scaler = c("silverman", "IQR", "sd", "var")){
  if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")
  
  nstar <- processnstar(nstar,pp)
  scaler <- processscaler(scaler,pp)
  
  RK <- 1/(4*pi)
  d <- 2
  V <- (16*gamma((d+8)/2)*d*(d+2))/((d+8)^((d+6)/2)*pi^(d/2))
  
  return(scaler*(((RK*d)/(nstar*V))^(1/(d+4))))
}