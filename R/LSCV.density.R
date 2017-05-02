LSCV.density <- function(pp,hlim=NULL,hseq=NULL,resolution=64,edge=TRUE,auto.optim=TRUE,seqres=30,parallelise=NULL,verbose=TRUE,type="spatial",lambdalim=NULL,lambdaseq=NULL,tlim=NULL){
  if(!is.null(hlim)){
    if(hlim[1]>=hlim[2]) stop("invalid h limits")
  }
  if(!is.null(lambdalim)){
    if(lambdalim[1]>=lambdalim[2]) stop("invalid lambda limits")
  }
  if(class(pp)!="ppp") stop("data object 'pp' must be of class \"ppp\"")
  W <- Window(pp)
  
  if(is.null(hlim)){
    md <- min(nndist(unique(pp)))
    hlim <- c(md,max(md*50,min(diff(W$xrange),diff(W$yrange))/6))
  }
  
  #if(type=="spatial"){
    if(auto.optim){
      if(verbose) cat("Searching for optimal h in [",round(hlim[1],3),",",round(hlim[2],3),"]...",sep="")
      result <- optimise(LSCV.density.spatial.single,interval=hlim,pp=pp,res=resolution,edge=edge)$minimum
      if(verbose) cat("Done.\n")
    } else {
      if(is.null(hseq)) hseq <- seq(hlim[1],hlim[2],length=seqres)
      hn <- length(hseq)
      if(is.null(parallelise)){
        lscv.vec <- rep(NA,hn)
        if(verbose) pb <- txtProgressBar(1,hn)
        for(i in 1:hn){
          lscv.vec[i] <- LSCV.density.spatial.single(hseq[i],pp,resolution,edge)
          if(verbose) setTxtProgressBar(pb,i)
        }
        if(verbose) close(pb)
      } else {
        ncores <- detectCores()
        if(verbose) cat(paste("Evaluating criterion on",parallelise,"/",ncores,"cores..."))
        if(parallelise>ncores) stop("cores requested exceeds available count")
        registerDoParallel(cores=parallelise)
        lscv.vec <- foreach(i=1:hn,.packages="spatstat",.combine=c) %dopar% {
          return(LSCV.density.spatial.single(hseq[i],pp,resolution,edge))
        }
        if(verbose) cat("Done.\n")
      }
      result <- cbind(hseq,lscv.vec)
      dimnames(result)[[2]] <- c("h","CV")
    }
  
  return(result)
  
  # } else if(type=="spattemp"){
  #   stop("'type' must be provided as \"spatial\"; all else currently unimplemented")
  #   if(auto.optim){
  #     strt <- c(NS(pp),bw.nrd0(marks(pp)))
  #     return(optim(par=strt,LSCV.density.spattemp.single,pp=pp,res=res,tlim=tlim,edge=edge,...)$par)
  #   } else {
  #     if(is.null(hseq)) hseq <- seq(hlim[1],hlim[2],length=seqres)
  #     hn <- length(hseq)
  #     if(is.null(lambdaseq)) lambdaseq <- seq(lambdalim[1],lambdalim[2],length=seqres)
  #     ln <- length(lambdaseq)
  #     hl <- expand.grid(hseq,lambdaseq)
  #     if(is.na(parallelise)){
  #       lscv.vec <- rep(NA,hn*ln)
  #       for(i in 1:(hn*ln)) lscv.vec[i] <- LSCV.density.spattemp.single(as.numeric(hl[i,]),pp,res,tlim,edge)
  #     } else {
  #       if(parallelise>detectCores()) stop("cores requested exceeds available count")
  #       registerDoParallel(cores=parallelise)
  #       lscv.vec <- foreach(i=1:(hn*ln),.packages="spatstat",.combine=c) %dopar% {
  #         return(LSCV.density.spattemp.single(as.numeric(hl[i,]),pp,res,tlim,edge))
  #       }
  #     }
  #     return(list(hs=hl[,1],ls=hl[,2],CV=lscv.vec))
  #   }
  # } else {
  #   stop("'type' must be provided as either \"spatial\" (spatial only) or \"spattemp\" (spatiotemporal)")
  # }
}



