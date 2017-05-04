
tol.mc.ada <- function(rs,ITER,parallel,verbose,...){
  fd <- rs$f
  gd <- rs$g
  pool <- suppressWarnings(superimpose(fd$pp,gd$pp))
  nf <- npoints(fd$pp)
  ng <- npoints(gd$pp)
  res <- dim(fd$z)[1]
  indx <- 1:npoints(pool)
  
  if(is.null(fd$q)){
    edg <- "none"
  } else if(is.im(fd$q)){
    edg <- "uniform"
  } else {
    edg <- "diggle"
  }
  
  logt <- !all(rs$rr>=0)
  rmat <- as.matrix(rs$rr)
  mcmat <- matrix(1,res,res)
  
  if(is.null(parallel)){
    if(verbose) pb <- txtProgressBar(0,ITER-1,style=3)
    for(i in 1:(ITER-1)){
      shuff <- sample(indx)
      rtemp <- as.matrix(risk(pool[shuff[1:nf]],pool[shuff[(nf+1):(nf+ng)]],log=logt,verbose=FALSE,h0=c(fd$h0,gd$h0),hp=c(fd$hp,gd$hp),resolution=res,edg=edg,adapt=TRUE,...)$rr)
      mcmat <- mcmat+(rtemp>=rmat)
      if(verbose) setTxtProgressBar(pb,i)
    }
    if(verbose) close(pb)
  } else {
    ncores <- detectCores()
    if(verbose) cat(paste("Running MC iterations on",parallel,"/",ncores,"cores..."))
    if(parallel>ncores) stop("Parallel cores requested exceeds available count")
    registerDoParallel(cores=parallel)
    mclist <- foreach(i=1:(ITER-1),.packages=c("spatstat","sparr")) %dopar% {
      shuff <- sample(indx)
      rtemp <- as.matrix(risk(pool[shuff[1:nf]],pool[shuff[(nf+1):(nf+ng)]],log=logt,verbose=FALSE,h0=c(fd$h0,gd$h0),hp=c(fd$hp,gd$hp),resolution=res,edg=edg,adapt=TRUE,...)$rr)
      return(rtemp>=rmat)
    }
    
    if(verbose) cat("Done.\n")
    mcmat <- Reduce("+",mclist) + mcmat
  }
  return(mcmat/ITER)
}



