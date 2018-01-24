#' @rdname CV
#' @export
LIK.density <- function(pp,hlim=NULL,hseq=NULL,resolution=64,edge=TRUE,auto.optim=TRUE,
                         type=c("fixed","adaptive"),seqres=30,parallelise=NULL,
                         zero.action=0,verbose=TRUE,...){
  
  if(class(pp)!="ppp") stop("data object 'pp' must be of class \"ppp\"")
  W <- Window(pp)
  
  if(is.null(hlim)){
    ppu <- pp
    marks(ppu) <- NULL
    md <- min(nndist(unique(ppu)))
    hlim <- c(md,max(md*50,min(diff(W$xrange),diff(W$yrange))/6))
  } else {
    hlim <- checkran(hlim,"'hlim'")
  }
  
  if(!zero.action%in%((-1):2)) stop("invalid 'zero.action'")
  
  typ <- type[1]
  if(typ=="fixed"){
    if(auto.optim){
      if(verbose) cat("Searching for optimal h in ",prange(hlim),"...",sep="")
      result <- suppressWarnings(optimise(LIK.density.spatial.single,interval=hlim,pp=pp,res=resolution,edge=edge,za=zero.action,maximum=TRUE)$maximum)
      if(verbose) cat("Done.\n")
    } else {
      if(is.null(hseq)) hseq <- seq(hlim[1],hlim[2],length=seqres)
      hn <- length(hseq)
      if(is.null(parallelise)){
        lik.vec <- rep(NA,hn)
        if(verbose) pb <- txtProgressBar(1,hn)
        for(i in 1:hn){
          lik.vec[i] <- LIK.density.spatial.single(hseq[i],pp,resolution,edge,za=zero.action)
          if(verbose) setTxtProgressBar(pb,i)
        }
        if(verbose) close(pb)
      } else {
        ncores <- detectCores()
        if(verbose) cat(paste("Evaluating criterion on",parallelise,"/",ncores,"cores..."))
        if(parallelise>ncores) stop("cores requested exceeds available count")
        registerDoParallel(cores=parallelise)
        lik.vec <- foreach(i=1:hn,.packages="spatstat",.combine=c) %dopar% {
          return(LIK.density.spatial.single(hseq[i],pp,resolution,edge,zero.action))
        }
        if(verbose) cat("Done.\n")
      }
      result <- cbind(hseq,lik.vec)
      dimnames(result)[[2]] <- c("h","CV")
    }
  } else if(typ=="adaptive"){
    
    ellip <- list(...)
    
    if(is.null(ellip$hp)){
      if(verbose) cat("Selecting pilot bandwidth...")
      hp <- LSCV.density(pp,verbose=FALSE,zero.action=zero.action)
      if(verbose) cat(paste("Done.\n   [ Found hp =",hp,"]\n"))
    } else {
      hp <- ellip$hp
    }
    
    if(is.null(ellip$pilot.density)){
      pilot.density <- pp
    } else {
      pilot.density <- ellip$pilot.density
    }
    
    if(is.null(ellip$gamma.scale)){
      gamma.scale <- "geometric"
    } else {
      gamma.scale <- ellip$gamma.scale
    }
    
    if(is.null(ellip$trim)){
      trim <- 5
    } else {
      trim <- ellip$trim
    }
    
    if(is.null(ellip$dimz)){
      dimz <- resolution
    } else {
      dimz <- ellip$dimz
    }
    
    if(verbose) cat("Computing multi-scale estimate...")
    hhash <- mean(hlim)
    msobject <- multiscale.density(pp,h0=hhash,hp=hp,h0fac=hlim/hhash,edge=ifelse(edge,"uniform","none"),resolution=resolution,dimz=dimz,gamma.scale=gamma.scale,trim=trim,intensity=FALSE,pilot.density=pilot.density,verbose=FALSE)
    if(verbose) cat("Done.\n")
    
    h0range <- range(as.numeric(names(msobject$z)))
    if(auto.optim){
      if(verbose) cat("Searching for optimal h0 in ",prange(h0range),"...",sep="")
      h0opt <- suppressWarnings(optimise(ms.loo.lik,interval=h0range,object=msobject,za=zero.action,maximum=TRUE)$maximum)
      if(verbose) cat("Done.\n")
      return(h0opt)
    } else {
      if(is.null(hseq)) hseq <- seq(h0range[1],h0range[2],length=seqres)
      hn <- length(hseq)
      lik.vec <- rep(NA,hn)
      if(verbose) pb <- txtProgressBar(1,hn)
      for(i in 1:hn){
        lik.vec[i] <- ms.loo.lik(hseq[i],msobject,zero.action)
        if(verbose) setTxtProgressBar(pb,i)
      }
      if(verbose) close(pb)
      
      result <- cbind(hseq,lik.vec)
      dimnames(result)[[2]] <- c("h0","CV")
    }
  } else stop("invalid 'type'")
  
  return(result)
}



