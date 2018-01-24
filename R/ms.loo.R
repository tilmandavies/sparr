ms.loo <- function(h0,object,za){
  X <- object$pp
  n <- npoints(X)
  
  requested <- multiscale.slice(object,h0,checkargs=FALSE)
  rh <- requested$h
  rz <- requested$z
  rq <- requested$q
  
  if(za==-1){
    if(any(rz<=0)) return(Inf)
  }
  
  rint <- integral(requested$z)
  zpoints <- safelookup(rz,X,warn=FALSE)
  
  if(is.null(rq)) qpoints <- rep(1,n)
  else qpoints <- safelookup(rq,X,warn=FALSE)
  
  rzn <- (rz/rint)^2
  loo.atpoints <- (zpoints-dnorm(0,sd=rh)^2/qpoints)/(n-1)

  rznint <- integral(rzn)
  if(any(loo.atpoints<=0)){
    if(za==2){
      loo.atpoints[loo.atpoints<=0] <- min(loo.atpoints[loo.atpoints>0])
    } else if(za==1){
      loo.atpoints <- posifybivden(loo.atpoints)
    } else {
      return(Inf)
    }
  } #was: return(rznint)
  return(rznint-2*mean(loo.atpoints))
}

ms.loo.lik <- function(h0,object,za){
  X <- object$pp
  n <- npoints(X)

  requested <- multiscale.slice(object,h0,checkargs=FALSE)
  rh <- requested$h
  rz <- requested$z
  rq <- requested$q
  
  if(za==-1){
    if(any(rz<=0)) return(-Inf)
  }
  
  zpoints <- safelookup(rz,X,warn=FALSE)
    
  if(is.null(rq)) qpoints <- rep(1,n)
  else qpoints <- safelookup(rq,X,warn=FALSE)
    
  loo.atpoints <- zpoints-(1/n)*(dnorm(0,sd=rh)^2/qpoints)
    
  if(any(loo.atpoints<=0)){
    if(za==2){
      loo.atpoints[loo.atpoints<=0] <- min(loo.atpoints[loo.atpoints>0])
    } else if(za==1){
      loo.atpoints <- posifybivden(loo.atpoints)
    } else {
      return(-Inf)
    }
  } #was: return(log(min(loo.atpoints[loo.atpoints>0])))
  
    
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
  # if(any(frz<=0)||any(grz<=0)) return(Inf) # pretty strict protection; generates optimise warnings by default
  
  frz[frz<=0] <- limz
  grz[grz<=0] <- limz

  f.fpoints <- safelookup(frz,fX,warn=FALSE)
  g.gpoints <- safelookup(grz,gX,warn=FALSE)
  f.gpoints <- safelookup(frz,gX,warn=FALSE)
  g.fpoints <- safelookup(grz,fX,warn=FALSE)

  if(is.null(frq)){
    fqpoints <- rep(1,n1)
    gqpoints <- rep(1,n2)
  } else {
    fqpoints <- safelookup(frq,fX,warn=FALSE)
    gqpoints <- safelookup(grq,gX,warn=FALSE)
  }
  
  fminus <- f.fpoints - dnorm(0,sd=frh)^2/n1/fqpoints
  fminus[fminus<=0] <- limz
  gminus <- g.gpoints - dnorm(0,sd=grh)^2/n2/gqpoints
  gminus[gminus<=0] <- limz
  
  if(!hazey) return(2*mean((log(f.gpoints) - log(gminus))/gminus) - 2*mean((log(fminus) - log(g.fpoints))/fminus) - integral((log(frz)-log(grz))^2))
  else return(mean((f.gpoints/gminus)^2) - 2*mean(fminus/g.fpoints))
}
