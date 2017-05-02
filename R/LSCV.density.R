#' Leave-one-out least-squares cross-validation (LSCV) bandwidth selector
#' 
#' Isotropic fixed bandwidth selection for standalone 2D density/intensity
#' based on classical unbiased cross-validation
#' 
#' 
#' @param pp An object of class \code{\link[spatstat]{ppp}} giving the observed
#'   2D data to be smoothed.
#' @param hlim An optional vector of length 2 giving the limits of the
#'   optimisation routine with respect to the bandwidth. If unspecified, the
#'   function attempts to choose this automatically.
#' @param hseq An optional increasing sequence of bandwidth values at which to
#'   manually evaluate the optimisation criterion. Used only in the case
#'   \code{(!auto.optim && is.null(hlim))}.
#' @param resolution Spatial grid size; the optimisation will be based on a
#'   [\code{resolution} \eqn{\times}{x} \code{resolution}] density estimate.
#' @param edge Logical value indicating whether to edge-correct the density
#'   estimates used.
#' @param auto.optim Logical value indicating whether to automate the numerical
#'   optimisation using \code{\link{optimise}}. If \code{FALSE}, the optimisation
#'   criterion is evaluated over \code{hseq} (if supplied), or over a seqence of
#'   values controlled by \code{hlim} and \code{seqres}.
#' @param seqres Optional resolution of an increasing sequence of bandwidth
#'   values. Only used if \code{(!auto.optim && is.null(hseq))}.
#' @param parallelise Numeric argument to invoke parallel processing, giving
#'   the number of CPU cores to use when \code{!auto.optim}. Experimental. Test
#'   your system first using \code{parallel::detectCores()} to identify the
#'   number of cores available to you.
#' @param verbose Logical value indicating whether to provide function progress
#'   commentary.
#' @param type Unimplemented. Future argument for spatiotemporal bandwidth
#'   selection.
#' @param lambdalim Unimplemented. Future argument for spatiotemporal bandwidth
#'   selection.
#' @param lambdaseq Unimplemented. Future argument for spatiotemporal bandwidth
#'   selection.
#' @param tlim Unimplemented. Future argument for spatiotemporal bandwidth
#'   selection.
#'
#' @return A single numeric value of the estimated bandwidth (if
#'   \code{auto.optim = TRUE}). Otherwise, a list of two numeric vectors of equal
#'   length giving the bandwidth sequence (as \code{hs}) and corresponding CV
#'   function value (as \code{CV}).
#'
#' @section Warning: Leave-one-out LSCV for bandwidth selection in kernel
#' density estimation is notoriously unstable in practice and has a tendency to
#' produce rather small bandwidths. Satisfactory bandwidths are not guaranteed
#' for every application. This method can also be computationally expensive for
#' large data sets and fine evaluation grid resolutions. The user may need to
#' experiment with adjusting \code{hlim} to find a suitable minimum.
#'
#' @author T. M. Davies
#'
#' @seealso Functions for bandwidth selection in package
#'   \code{\link{spatstat}}: \code{\link[spatstat]{bw.diggle}};
#'   \code{\link[spatstat]{bw.ppl}}; \code{\link[spatstat]{bw.scott}};
#'   \code{\link[spatstat]{bw.frac}}.
#'
#' @references
#' Silverman, B.W. (1986), \emph{Density Estimation for Statistics
#' and Data Analysis}, Chapman & Hall, New York.
#'
#' Wand, M.P. and Jones,
#' C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.
#'
#' @examples
#' 
#' ## To be filled
#' 
LSCV.density <- function(pp,hlim=NULL,hseq=NULL,resolution=64,edge=TRUE,auto.optim=TRUE,
                         seqres=30,parallelise=NULL,verbose=TRUE,type="spatial",
                         lambdalim=NULL,lambdaseq=NULL,tlim=NULL){
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



