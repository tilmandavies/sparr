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
  loo.atpoints <- (zpoints-dnorm(0,sd=rh)^2/qpoints)/(n-1)

  rznint <- integral(rzn)
  if(any(loo.atpoints<=0)) return(rznint)
  return(rznint-2*mean(loo.atpoints))
}

ms.loo.lik <- function(h0,object){
  X <- object$pp
  n <- npoints(X)

  requested <- multiscale.slice(object,h0,checkargs=FALSE)
  rh <- requested$h
  rz <- requested$z
  rq <- requested$q
  zpoints <- safelookup(rz,X,warn=FALSE)
    
  if(is.null(rq)) qpoints <- rep(1,n)
  else qpoints <- safelookup(rq,X,warn=FALSE)
    
  loo.atpoints <- zpoints-(1/n)*(dnorm(0,sd=rh)^2/qpoints)
    
  if(any(loo.atpoints<=0)) return(log(min(loo.atpoints[loo.atpoints>0])))
  return(mean(log(loo.atpoints)))
} 

ms.loo.risk <- function(h0,fob,gob,hazey=FALSE){
  fX <- fob$pp
  gX <- gob$pp
  n1 <- npoints(fX)
  n2 <- npoints(gX)
  
  f.requested <- multiscale.slice(fob,h0,checkargs=FALSE)
  g.requested <- multiscale.slice(gob,h0,checkargs=FALSE)
  
  frh <- f.requested$h
  grh <- g.requested$h
  frz <- f.requested$z
  grz <- g.requested$z
  frq <- f.requested$q
  grq <- g.requested$q
  
  limz <- min(c(min(frz[frz>0]),min(grz[grz>0])))
  f.fpoints <- safelookup(frz,fX,warn=FALSE)
  f.fpoints[f.fpoints<=0] <- limz
  g.gpoints <- safelookup(grz,gX,warn=FALSE)
  g.gpoints[g.gpoints<=0] <- limz
  f.gpoints <- safelookup(frz,gX,warn=FALSE)
  f.gpoints[f.gpoints<=0] <- limz
  g.fpoints <- safelookup(grz,fX,warn=FALSE)
  g.fpoints[g.fpoints<=0] <- limz

  if(is.null(frq)){
    fqpoints <- rep(1,n1)
    gqpoints <- rep(1,n2)
  } else {
    fqpoints <- safelookup(frq,fX,warn=FALSE)
    gqpoints <- safelookup(grq,gX,warn=FALSE)
  }
  
  fminus <- f.fpoints - dnorm(0,sd=frh)^2/n1/fqpoints
  fminus[fminus<=0] <- limz # small bw protector
  gminus <- g.gpoints - dnorm(0,sd=grh)^2/n2/gqpoints
  gminus[gminus<=0] <- limz # small bw protector
  
  if(!hazey) return(2*mean((log(f.gpoints) - log(gminus))/gminus) - 2*mean((log(fminus) - log(g.fpoints))/fminus) - integral((log(frz)-log(grz))^2))
  else return(mean((f.gpoints/gminus)^2) - 2*mean(fminus/g.fpoints))
}
