bivden.LOO <- function(pp,h0,hp,edge,gamma.scale,trim,resolution,parallel,weights,za){
  n <- npoints(pp)
  if(is.null(weights)) weights <- rep(1,n)
  pilot.density.spec.loo <- density(pp,sigma=hp,dimyx=rep(resolution,2),at="points",edge=edge,positive=TRUE,leaveoneout=TRUE,weights=weights/(n-1))
  
  W <- Window(pp)
  Wm <- as.mask(W,dimyx=rep(resolution,2))
  evalxy <- as.matrix(expand.grid(Wm$xcol,Wm$yrow))
  notin <- !inside.owin(x=evalxy[,1],y=evalxy[,2],w=W)
  evalxy.in <- evalxy[!notin,]
  
  hsi <- 1:n
  loo <- rep(NA,n)
  qv <- rep(1,n)
  h.spec <- matrix(NA,n,n)
  if(is.null(parallel)){
    for(i in hsi){
      ppmi <- pp[-i]
      wi <- weights[-i]
      
      pilot.density <- density(ppmi,sigma=hp,dimyx=rep(resolution,2),edge=edge,positive=FALSE,weights=wi)
      pilot.density.spec <- safelookup(pilot.density,ppmi,warn=FALSE)#density(ppmi,sigma=hp,dimyx=rep(resolution,2),at="points",positive=TRUE,leaveoneout=FALSE)
      pi.int <- integral(pilot.density)
      pilot.density <- pilot.density/pi.int
      pilot.density.spec <- pilot.density.spec/pi.int
      
      if(za==0){
        if(any(pilot.density.spec<=0)){
          loo[i] <- NA
          next
        }
      } else if(za==1){
        pilot.density.spec <- posifybivden(pilot.density.spec)
      } else if(za==2){
        pilot.density.spec[pilot.density.spec<=0] <- min(pilot.density.spec[pilot.density.spec>0])
      }
      
      pspec <- pilot.density.spec^(-0.5)
      gamma <- processgamma(gamma.scale,pilot.density.spec)
      
      h.spec.mi <- h.spec[hsi[-i],i] <- h0*pmin(pspec/gamma,trim)
      h.hypo.i <- h.spec[i,i] <- h0*min(pilot.density.spec.loo[i]^(-0.5)/gamma,trim)
      
      if(edge){
        gxy <- kernel2d(evalxy.in[,1]-pp$x[i],evalxy.in[,2]-pp$y[i],h.hypo.i)
        qv[i] <- dintegral(gxy,Wm$xstep,Wm$ystep)
      }
      ivals <- kernel2d(pp$x[i]-ppmi$x, pp$y[i]-ppmi$y,h.spec.mi)
      loo[i] <- mean(wi*ivals)/qv[i]  #min(,infvec[i])
    }
  } else {
    if(parallel>detectCores()) stop("Parallel cores requested exceeds available count")
    
    registerDoParallel(cores=parallel)
    loo <- foreach(i=1:n,.packages="spatstat",.combine=c) %dopar% {
      ppmi <- pp[-i]
      wi <- weights[-i]
      
      pilot.density <- density(ppmi,sigma=hp,dimyx=rep(resolution,2),edge=edge,positive=TRUE,weights=wi)
      pilot.density.spec <- safelookup(pilot.density,ppmi,warn=FALSE)#density(ppmi,sigma=hp,dimyx=rep(resolution,2),at="points",positive=TRUE,leaveoneout=FALSE)
      pi.int <- integral(pilot.density)
      pilot.density$v <- pilot.density$v/pi.int
      pilot.density.spec <- pilot.density.spec/pi.int
      
      if(za==0){
        if(any(pilot.density.spec<=0)){
          return(NA)
        }
      } else if(za==1){
        pilot.density.spec <- posifybivden(pilot.density.spec)
      } else if(za==2){
        pilot.density.spec[pilot.density.spec<=0] <- min(pilot.density.spec[pilot.density.spec>0])
      }
      
      pspec <- pilot.density.spec^(-0.5)
      gamma <- processgamma(gamma.scale,pilot.density.spec)
      
      h.spec.mi <- h0*pmin(pspec/gamma,trim) # h.spec[hsi[-i],i] <- 
      h.hypo.i <- h0*min(pilot.density.spec.loo[i]^(-0.5)/gamma,trim) # h.spec[i,i] <- 
      
      if(edge){
        gxy <- kernel2d(evalxy.in[,1]-pp$x[i],evalxy.in[,2]-pp$y[i],h.hypo.i)
        qi <- dintegral(gxy,Wm$xstep,Wm$ystep)
      } else {
        qi <- 1
      }
      ivals <- kernel2d(pp$x[i]-ppmi$x,pp$y[i]-ppmi$y,h.spec.mi)
      
      return(mean(wi*ivals)/qi)
    }	
  }
  return(list(loo,qv,h.spec)) # qv and h.spec only filled in when parallel=NULL
}