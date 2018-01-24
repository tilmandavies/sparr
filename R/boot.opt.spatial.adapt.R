boot.opt.spatial.adapt <- function(pp,h0ref,h0fac,hp,edg,refden,N,B,res,dimz,verbose,parallelise,...){
    
  if(is.na(parallelise)){
    isemat <- matrix(NA,N,B)
    if(verbose) pb <- txtProgressBar(0,N)
    for(i in 1:N){
      tempdata <- rimpoly(pp$n,refden,Window(pp))
      tempadapt <- multiscale.density(tempdata,h0=h0ref,hp=hp,h0fac=h0fac,edge=edg,resolution=res,dimz=dimz,verbose=FALSE,...)
      h0seq <- seq(tempadapt$h0range[1],tempadapt$h0range[2],length=B)
      for(j in 1:B){
        bj <- multiscale.slice(tempadapt,h0seq[j])
        isemat[i,j] <- integral((bj$z-refden)^2)
      }
      if(verbose) setTxtProgressBar(pb,i)
    }
    if(verbose) close(pb)
    resultmat <- cbind(h0seq,colMeans(isemat))
    
  } else {
    
    totcor <- detectCores()
    if(parallelise>totcor) stop("Parallel cores requested exceeds available count")
    
    if(verbose) cat(paste("   --bootstrapping on",parallelise,"/",totcor,"cores--\n"))
    registerDoParallel(cores=parallelise)
    isemat <- foreach(i=1:N,.packages=c("sparr","spatstat"),.combine=rbind) %dopar% {
      tempdata <- rimpoly(pp$n,refden,Window(pp))
      tempadapt <- multiscale.density(tempdata,h0=h0ref,hp=hp,h0fac=h0fac,edge=edg,resolution=res,dimz=dimz,verbose=FALSE,...)
      h0seq <- seq(tempadapt$h0range[1],tempadapt$h0range[2],length=B)
      isevec <- rep(NA,B)
      for(j in 1:B){
        bj <- multiscale.slice(tempadapt,h0seq[j])
        isevec[j] <- integral((bj$z-refden)^2)
      }
      return(rbind(isevec,h0seq))
    }
    if(verbose) cat("Done.\n")
    h0seq <- isemat[2,]
    isemat <- isemat[seq(1,2*N-1,2),]
    resultmat <- cbind(h0seq,colMeans(isemat))
  }
  
  #old return val: resultmat[which.min(resultmat[,2]),1][1]
  rs <- spline(resultmat)
  return(list(h=rs$x[which.min(rs$y)],mat=resultmat))
}