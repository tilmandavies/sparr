spattemp.LOO <- function(pp,tt,h,lambda,tlim,xyin,xys,sedge,tedge,parallelise){
  W <- Window(pp)
  n <- npoints(pp)
  ppmat <- cbind(pp$x,pp$y)
  qs <- qt <- 1
  
  if(is.na(parallelise)){
    loo <- rep(NA,n)
    for(i in 1:n){
      ppt.i <- pp[i]
      ppt.mi <- pp[-i]
      t.i <- tt[i]
      t.mi <- tt[-i]
      
      if(sedge){
        pxy <- kernel2d(xyin[,1]-ppt.i$x, xyin[,2]-ppt.i$y, h)
        qs <- dintegral(pxy,xys[1],xys[2])
      }
      if(tedge) qt <- pnorm(tlim[2],t.i,lambda) - pnorm(tlim[1],t.i,lambda)
      
      ut <- (t.i-t.mi)/lambda
      ivals <- kernel2d(ppt.i$x-ppt.mi$x,ppt.i$y-ppt.mi$y,h)*lambda^(-1)*exp(-0.5*ut^2)/sqrt(2*pi)
      loo[i] <- mean(ivals)/(qs*qt)
    }
  } else {
    registerDoParallel(cores=parallelise)
    loo <- foreach(i=1:n,.packages="spatstat",.combine=c) %dopar% {
      ppt.i <- pp[i]
      ppt.mi <- pp[-i]
      t.i <- tt[i]
      t.mi <- tt[-i]
      
      if(sedge){
        pxy <- kernel2d(xyin[,1]-ppt.i$x, xyin[,2]-ppt.i$y, h)
        qs <- dintegral(pxy,xys[1],xys[2])
      }
      if(tedge) qt <- pnorm(tlim[2],t.i,lambda) - pnorm(tlim[1],t.i,lambda)
      
      ut <- (t.i-t.mi)/lambda
      ivals <- kernel2d(ppt.i$x-ppt.mi$x,ppt.i$y-ppt.mi$y,h)*lambda^(-1)*exp(-0.5*ut^2)/sqrt(2*pi)
      return(mean(ivals)/(qs*qt))
    }
  }
  return(loo)
}