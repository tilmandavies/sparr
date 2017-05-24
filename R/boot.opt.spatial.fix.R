boot.opt.spatial.fix <- function(h,rmdiag,edg,WM,resolution,fM,ifft_scale,inn,GN,evalyx.redu,
                                 pp,epsilon.eta,eta,nn,boot3,use_fftw,parallelise){
  
  if(edg){
    # get edge surface for 'h'
    fK.h <- kernel2d_fft(h,WM$xstep,WM$ystep,resolution)
    fK.con <- fft2d(fM*fK.h,inverse=TRUE,use_fftw)[1:resolution,1:resolution]
    edg.h <- Mod(fK.con)*ifft_scale
    edg.h[edg.h>1] <- 1  
    GE.h <- edg.h[inn]
  } else {
    GE.h <- rep(1,length(inn))
  }
  
  # do edge-corrected fixed-bandwidth spatial bootstrap for 'h'; no resampling necessary
  if(!rmdiag){
    
    if(is.na(parallelise)){
      b1c <- b2c <- rep(NA,GN)
      for(i in 1:GN){
        evx <- evalyx.redu[i,2]-pp$x
        evy <- evalyx.redu[i,1]-pp$y
        b2temp <- sum(epsilon.eta^(-1)*kernel2d(evx,evy,sqrt(h^2+eta^2)))
        b2c[i] <- b2temp*sum(epsilon.eta^(-1)*kernel2d(evx,evy,eta))
        b1c[i] <- (2*sqrt(pi)*h)^(-2)*sum(epsilon.eta^(-1)*kernel2d(evx,evy,sqrt(0.5*h^2+eta^2))) + ((nn-1)/nn)*b2temp^2
      }
      boot1 <- GE.h^(-2)*b1c
      boot2 <- -2*GE.h^(-1)*b2c
    } else {
      totcor <- detectCores()
      if(parallelise>totcor) stop("Parallel cores requested exceeds available count")
      # if(verbose) cat(paste("   --optimising on",parallelise,"/",totcor,"cores--\n"))
      registerDoParallel(cores=parallelise)
      b12c <- foreach(i=1:GN,.packages="sparr",.combine=rbind) %dopar% {
        evx <- evalyx.redu[i,2]-pp$x
        evy <- evalyx.redu[i,1]-pp$y
        b2temp <- sum(epsilon.eta^(-1)*kernel2d(evx,evy,sqrt(h^2+eta^2)))
        b2c <- b2temp*sum(epsilon.eta^(-1)*kernel2d(evx,evy,eta))
        b1c <- (2*sqrt(pi)*h)^(-2)*sum(epsilon.eta^(-1)*kernel2d(evx,evy,sqrt(0.5*h^2+eta^2))) + ((nn-1)/nn)*b2temp^2
        result <- c(b1c,b2c)
      }
      boot1 <- GE.h^(-2)*b12c[,1]
      boot2 <- -2*GE.h^(-1)*b12c[,2]
    }
    
    return(sum((nn^(-2)*boot1+nn^(-2)*boot2+boot3)*WM$xstep*WM$ystep))
    
  } else {
    
    if(is.na(parallelise)){
      b1c <- b2c <- rep(NA,GN)
      for(i in 1:GN){
        evx <- evalyx.redu[i,2]-pp$x
        evy <- evalyx.redu[i,1]-pp$y
        b2temp.a <- epsilon.eta^(-1)*kernel2d(evx,evy,sqrt(h^2+eta^2))
        b2temp.b <- epsilon.eta^(-1)*kernel2d(evx,evy,eta)
        b2c[i] <- sum(b2temp.a)*sum(b2temp.b) - sum(b2temp.a*b2temp.b)
        b1c[i] <- (2*sqrt(pi)*h)^(-2)*sum(epsilon.eta^(-1)*kernel2d(evx,evy,sqrt(0.5*h^2+eta^2))) + ((nn-1)/nn)*sum(b2temp.a)^2
      }
      boot1 <- GE.h^(-2)*b1c
      boot2 <- -2*GE.h^(-1)*b2c
    } else {
      totcor <- detectCores()
      if(parallelise>totcor) stop("Parallel cores requested exceeds available count")
      # if(verbose) cat(paste("   --optimising on",parallelise,"/",totcor,"cores--\n"))
      registerDoParallel(cores=parallelise)
      b12c <- foreach(i=1:GN,.packages="sparr",.combine=rbind) %dopar% {
        evx <- evalyx.redu[i,2]-pp$x
        evy <- evalyx.redu[i,1]-pp$y
        b2temp.a <- epsilon.eta^(-1)*kernel2d(evx,evy,sqrt(h^2+eta^2))
        b2temp.b <- epsilon.eta^(-1)*kernel2d(evx,evy,eta)
        b2c <- sum(b2temp.a)*sum(b2temp.b) - sum(b2temp.a*b2temp.b)
        b1c <- (2*sqrt(pi)*h)^(-2)*sum(epsilon.eta^(-1)*kernel2d(evx,evy,sqrt(0.5*h^2+eta^2))) + ((nn-1)/nn)*sum(b2temp.a)^2
        return(c(b1c,b2c))
      }    
      boot1 <- GE.h^(-2)*b12c[,1]
      boot2 <- -2*GE.h^(-1)*b12c[,2]
    }
    
    return(sum(nn^(-2)*(boot1+boot2)*WM$xstep*WM$ystep)+boot3)
    
  }
}