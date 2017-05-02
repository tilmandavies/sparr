LSCV.risk <- function(f, g = NULL, hlim = NULL, hseq = NULL, method = c("kelsall-diggle", "hazelton", "davies"), resolution = 64, edge = TRUE, auto.optim = TRUE, seqres = 30, parallelise = NULL, verbose = TRUE){
  if(!inherits(f,"ppp")) stop("'f' must be an object of class \"ppp\"")
  if(is.null(g)){
    fm <- marks(f)
    if(!is.factor(fm)) marks(f) <- fm <- factor(fm)
    if(nlevels(fm)!=2) stop("'f' marks must be dichotomous if 'g' unsupplied")
    fs <- split(f)
    f <- fs[[1]]
    g <- fs[[2]]
  }
  if(!inherits(g,"ppp")) stop("'g' must be an object of class \"ppp\"")

  W <- Window(f)
  if(!identical_windows(W,Window(g))) stop("study windows for 'f' and 'g' must be identical")
  
  if(!is.null(hlim)){
    if(hlim[1]>=hlim[2]) stop("invalid 'hlim'")
  } else {
    md <- min(c(nndist(unique(f)),nndist(unique(g))))
    hlim <- c(md,max(md*50,min(diff(W$xrange),diff(W$yrange))/6))
  }

  meth <- method[1]
  
  if(auto.optim){
    if(meth=="kelsall-diggle"){
      if(verbose) cat("Searching for optimal Kelsall-Diggle h in [",round(hlim[1],3),",",round(hlim[2],3),"]...",sep="")
      result <- optimise(LSCV.risk.single,interval=hlim,cases=f,controls=g,res=resolution,edge=edge,hazey=FALSE)$minimum
    } else if(meth=="hazelton"){
      if(verbose) cat("Searching for optimal Hazelton h in [",round(hlim[1],3),",",round(hlim[2],3),"]...",sep="")
      result <- optimise(LSCV.risk.single,interval=hlim,cases=f,controls=g,res=resolution,edge=edge,hazey=TRUE)$minimum
    } else if(meth=="davies"){
      if(verbose) cat("Searching for optimal Davies h in [",round(hlim[1],3),",",round(hlim[2],3),"]\n  -initialisation...",sep="")
      pooled <- suppressWarnings(superimpose(f,g))
      lambda <- LSCV.density(pooled,verbose=FALSE)
      bp <- BAMprep(f,g,lambda,3,resolution)
      if(verbose) cat("Done.\n  -optimisation...")
      result <- optimise(BAM.single,interval=hlim,edge=edge,BP=bp)$minimum
    } else {
      stop("invalid 'method'")
    }
    if(verbose) cat("Done.\n")
  } else {
    if(is.null(hseq)) hseq <- seq(hlim[1],hlim[2],length=seqres)
    hn <- length(hseq)
    if(meth=="kelsall-diggle"){
      if(is.null(parallelise)){
        lscv.vec <- rep(NA,hn)
        if(verbose) pb <- txtProgressBar(1,hn)
        for(i in 1:hn){
          lscv.vec[i] <- LSCV.risk.single(hseq[i],cases=f,controls=g,res=resolution,edge=edge,hazey=FALSE)
          if(verbose) setTxtProgressBar(pb,i)
        }
        if(verbose) close(pb)
      } else {
        ncores <- detectCores()
        if(verbose) cat(paste("Evaluating criterion on",parallelise,"/",ncores,"cores..."))
        if(parallelise>ncores) stop("cores requested exceeds available count")
        registerDoParallel(cores=parallelise)
        lscv.vec <- foreach(i=1:hn,.packages="spatstat",.combine=c) %dopar% {
          return(LSCV.risk.single(hseq[i],cases=f,controls=g,res=resolution,edge=edge,hazey=FALSE))
        }
        if(verbose) cat("Done.\n")
      }
    } else if(meth=="hazelton"){
      if(is.null(parallelise)){
        lscv.vec <- rep(NA,hn)
        if(verbose) pb <- txtProgressBar(1,hn)
        for(i in 1:hn){
          lscv.vec[i] <- LSCV.risk.single(hseq[i],cases=f,controls=g,res=resolution,edge=edge,hazey=TRUE)
          if(verbose) setTxtProgressBar(pb,i)
        }
        if(verbose) close(pb)
      } else {
        ncores <- detectCores()        
        if(verbose) cat(paste("Evaluating criterion on",parallelise,"/",ncores,"cores..."))
        if(parallelise>ncores) stop("cores requested exceeds available count")
        registerDoParallel(cores=parallelise)
        lscv.vec <- foreach(i=1:hn,.packages="spatstat",.combine=c) %dopar% {
          return(LSCV.risk.single(hseq[i],cases=f,controls=g,res=resolution,edge=edge,hazey=TRUE))
        }
        if(verbose) cat("Done.\n")
      }
    } else if(meth=="davies"){
      pooled <- suppressWarnings(superimpose(f,g))
      lambda <- LSCV.density(pooled,verbose=FALSE)
      bp <- BAMprep(f,g,lambda,3,resolution)
      if(is.null(parallelise)){
        lscv.vec <- rep(NA,hn)
        if(verbose) pb <- txtProgressBar(1,hn)
        for(i in 1:hn){
          lscv.vec[i] <- BAM.single(hseq[i],edge=edge,BP=bp)
          if(verbose) setTxtProgressBar(pb,i)
        }
        if(verbose) close(pb)
      } else {
        ncores <- detectCores()
        if(verbose) cat(paste("Evaluating criterion on",parallelise,"/",ncores,"cores..."))
        if(parallelise>ncores) stop("cores requested exceeds available count")
        registerDoParallel(cores=parallelise)
        lscv.vec <- foreach(i=1:hn,.packages="spatstat",.combine=c) %dopar% {
          return(BAM.single(hseq[i],edge=edge,BP=bp))
        }
        if(verbose) cat("Done.\n")
      }
    } else {
      stop("invalid 'method'")
    }
    result <- cbind(hseq,lscv.vec)
    dimnames(result)[[2]] <- c("h","CV")
  }
  return(result)
}