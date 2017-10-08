tol.asy.st <- function(f,g,pooled,verbose){
  flen <- length(f$z)
  
  RK <- 1/(4*pi)
  RL <- rep(1/(2*sqrt(pi)),flen)
  
  nf <- npoints(f$pp)
  ng <- npoints(g$pp)

  if(verbose) message("   --convolution 1--")
  M <- Window(f$z[[1]])
  inside <- M$m
  pres <- nrow(inside)
  res2 <- 2*pres
  resseq <- 1:pres
  xcol.ker <- M$xstep*c(0:(pres-1),-rev(resseq))
  yrow.ker <- M$ystep*c(0:(pres-1),-rev(resseq))
  kerpixarea <- M$xstep*M$ystep
  len.pad <- res2^2
  Mpad <- matrix(0, ncol=2*pres, nrow=2*pres)
  Mpad[1:pres, 1:pres] <- inside
  fM <- fft(Mpad)
  
  qb <- function(oo,hfac,h){
    hfp <- hfac*h
    densX.ker <- dnorm(xcol.ker,sd=hfp)
    densY.ker <- dnorm(yrow.ker,sd=hfp)
    Kern <- outer(densY.ker,densX.ker,"*")*kerpixarea
    con <- fft(fM*fft(Kern), inverse=TRUE)/len.pad
    qhz <- im(Mod(con[1:pres,1:pres]),xcol=oo$spatial.z$xcol,yrow=oo$spatial.z$yrow)
    qhz[qhz>1] <- 1
    qhz[!inside] <- NA
    return(as.matrix(qhz))
  }
  
  if(verbose) message("   --convolution 2--")
  if(!is.null(pooled)){
    h <- pooled$h
    lam <- pooled$lambda
    
    qs <- as.matrix(pooled$qs)
    if(!all(qs==1)) RK <- RK*qb(pooled,sqrt(0.5),h)/qs^2
    
    qt <- pooled$qt
    nearedge <- which(qt<1)
    if(length(nearedge)>0){
      qt2 <- rep(1,length(qt))
      tgr <- as.numeric(names(pooled$z))
      for(i in nearedge) qt2[i] <- pnorm(pooled$tlim[2],mean=tgr[i],sd=lam/sqrt(2)) - pnorm(pooled$tlim[1],mean=tgr[i],sd=lam/sqrt(2))
      RL <- RL*qt2/qt^2
    }
    
    sig2 <- sig2.cond <- list()
    for(i in 1:flen){
      prefix <- (RK*RL[i])/(h^2*lam)
      sig2[[i]] <- prefix*(1/nf+1/ng)/as.matrix(pooled$z[[i]])
      sig2.cond[[i]] <- prefix*(1/(nf*f$temporal.z$y[i])+1/(ng*g$temporal.z$y[i]))/as.matrix(pooled$z.cond[[i]])
    }

  } else {
    h <- sqrt(prod(c(f$h,g$h0)))
    lam <- f$lambda
    qs <- as.matrix(f$qs)
    
    if(h!=f$h) qs <- qb(f,1,h)

    if(!all(qs==1)) RK <- RK*qb(f,sqrt(0.5),h)/qs^2
    
    qt <- f$qt
    nearedge <- which(qt<1)
    if(length(nearedge)>0){
      qt2 <- rep(1,length(qt))
      tgr <- as.numeric(names(f$z))
      for(i in nearedge) qt2[i] <- pnorm(f$tlim[2],mean=tgr[i],sd=lam/sqrt(2)) - pnorm(f$tlim[1],mean=tgr[i],sd=lam/sqrt(2))
      RL <- RL*qt2/qt^2
    }
    
    sig2 <- list()
    gadd <- RK/(as.matrix(g$z)*ng*h^2)
    for(i in 1:length(f$z)) sig2[[i]] <- (RK*RL[i])/(h^2*lam*nf*as.matrix(f$z[[i]])) + gadd
    sig2.cond <- sig2
  }
  
  
  # print(summary(as.vector(sig2[[40]])))
  
  
  return(list(v=sig2,vc=sig2.cond))
}
  