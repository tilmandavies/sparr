bivden.LOO <- function(pp,h0,hp,gamma.scale,trim,resolution,parallel){
  n <- npoints(pp)
  loo <- rep(NA,n)
  pilot.density.spec.loo <- density(pp,sigma=hp,dimyx=rep(resolution,2),at="points",positive=TRUE,leaveoneout=TRUE,weights=rep(1/(n-1),n))
  Wm <- as.mask(pp$window,dimyx=rep(resolution,2))
  
  if(is.null(parallel)){
    for(i in 1:n){
      ppmi <- pp[-i]
      
      pilot.density <- density(ppmi,sigma=hp,dimyx=rep(resolution,2),positive=TRUE)
      pilot.density.spec <- safelookup(pilot.density,ppmi,warn=FALSE)#density(ppmi,sigma=hp,dimyx=rep(resolution,2),at="points",positive=TRUE,leaveoneout=FALSE)
      pi.int <- integral(pilot.density)
      pilot.density <- pilot.density/pi.int
      pilot.density.spec <- pilot.density.spec/pi.int
      pspec <- pilot.density.spec^(-0.5)
      gamma <- processgamma(gamma.scale,pilot.density.spec)
      
      # PREVIOUS TRIMMING REGIMEN #
      # 			h.spec.mi <- h0*pilot.density.spec^(-0.5)/gamma
      #       beta.h <- trim*median(h.spec.mi,na.rm=TRUE)
      #       h.spec.mi[h.spec.mi>beta.h] <- beta.h
      #       h.hypo.i <- h0*pilot.density.spec.loo[i]^(-0.5)/gamma
      #       h.hypo.i[h.hypo.i>beta.h] <- beta.h
      
      # NEW TRIMMING REGIMEN #
      h.spec.mi <- h0*pmin(pspec/gamma,trim)
      #h.hypo.i <- h0*im(matrix(pmin(as.vector(as.matrix(pilot.density^(-0.5)))/gamma,trim),resolution,resolution),xcol=pilot.density$xcol,yrow=pilot.density$yrow)
      h.hypo.i <- h0*min(pilot.density.spec.loo[i]^(-0.5)/gamma,trim)
      
      uxy <- cbind(pp$x[i]-ppmi$x,pp$y[i]-ppmi$y)/h.spec.mi
      pr <- rnorm(1000*2,c(pp$x[i],pp$y[i]),h.hypo.i)
      indexer <- seq(1,1000*2,2)
      qhz <- mean(inside.owin(pr[indexer],pr[indexer+1],w=pp$window))
      
      loo[i] <- mean(h.spec.mi^(-2)*(exp(-0.5*rowSums(uxy^2))/(2*pi)))/qhz
    }
  } else {
    if(parallel>detectCores()) stop("Parallel cores requested exceeds available count")
    
    registerDoParallel(cores=parallel)
    loo <- foreach(i=1:n,.packages="spatstat",.combine=c) %dopar% {
      ppmi <- pp[-i]
      
      pilot.density <- density(ppmi,sigma=hp,dimyx=rep(resolution,2),positive=TRUE)
      pilot.density.spec <- safelookup(pilot.density,ppmi,warn=FALSE)#density(ppmi,sigma=hp,dimyx=rep(resolution,2),at="points",positive=TRUE,leaveoneout=FALSE)
      pi.int <- integral(pilot.density)
      pilot.density$v <- pilot.density$v/pi.int
      pilot.density.spec <- pilot.density.spec/pi.int
      pspec <- pilot.density.spec^(-0.5)
      gamma <- processgamma(gamma.scale,pilot.density.spec)
      
      h.spec.mi <- h0*pmin(pspec/gamma,trim)
      #h.hypo.i <- h0*im(matrix(pmin(as.vector(as.matrix(pilot.density^(-0.5)))/gamma,trim),resolution,resolution),xcol=pilot.density$xcol,yrow=pilot.density$yrow)
      h.hypo.i <- h0*min(pilot.density.spec.loo[i]^(-0.5)/gamma,trim)
      
      uxy <- cbind(pp$x[i]-ppmi$x,pp$y[i]-ppmi$y)/h.spec.mi
      pr <- rnorm(1000*2,c(pp$x[i],pp$y[i]),h.hypo.i)
      indexer <- seq(1,1000*2,2)
      qhz <- mean(inside.owin(pr[indexer],pr[indexer+1],w=pp$window))
      
      return(mean(h.spec.mi^(-2)*(exp(-0.5*rowSums(uxy^2))/(2*pi)))/qhz)
    }	
  }
  return(loo)
}