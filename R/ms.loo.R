ms.loo <- function(h0,object){
  X <- object$pp
  n <- npoints(X)
  
  requested <- multiscale.slice(object,h0,checkargs=FALSE)
  rh <- requested$h
  rz <- requested$z
  rq <- requested$q
  rint <- integral(requested$z)
  zpoints <- safelookup(rz,X,warn=FALSE)
  
  if(is.null(rq)) qpoints <- rep(1,n)
  else qpoints <- safelookup(rq,X,warn=FALSE)
  
  rzn <- (rz/rint)^2
  # rzn <- posifybivden(rzn)
  loo.atpoints <- (zpoints-(1/n)*dnorm(0,sd=rh)^2/qpoints)/(n-1)
  # loo.atpoints <- posifybivden(loo.atpoints)
  
  rznint <- integral(rzn)
  if(any(loo.atpoints<=0)) return(rznint)
  return(rznint-2*mean(loo.atpoints))
}