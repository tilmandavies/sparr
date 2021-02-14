#' Cross-validation bandwidths for spatial kernel density estimates
#' 
#' Isotropic fixed or global (for adaptive) bandwidth selection for standalone 2D density/intensity
#' based on either unbiased least squares cross-validation (LSCV) or likelihood (LIK) cross-validation.
#' 
#' This function implements the bivariate, edge-corrected versions of fixed-bandwidth least squares cross-validation and likelihood cross-validation
#' as outlined in Sections 3.4.3 and 3.4.4 of Silverman (1986) in order to select an optimal fixed smoothing bandwidth. With \code{type = "adaptive"} it may also be used to select the global bandwidth
#' for adaptive kernel density estimates, making use of multi-scale estimation (Davies and Baddeley, 2018) via \code{\link{multiscale.density}}.
#' Note that for computational reasons, the leave-one-out procedure is not performed on the pilot density in the adaptive setting; it
#' is only performed on the final stage estimate. Current development efforts include extending this functionality, see \code{\link{SLIK.adapt}}. See also `Warning' below.
#' 
#' Where \code{LSCV.density} is based on minimisation of an unbiased estimate of the mean integrated squared error (MISE) of the density, \code{LIK.density} is based on
#' maximisation of the cross-validated leave-one-out average of the log-likelihood of the density estimate with respect to \eqn{h}.
#' 
#' In both functions, the argument \code{zero.action} can be used to control the level of severity in response to small bandwidths that result (due to numerical error) in at least one density value being zero or less.
#' When \code{zero.action = -1}, the function strictly forbids bandwidths that would result in one or more \emph{pixel} values of a kernel estimate of the original data (i.e. anything over the whole region) being zero or less---this is the most restrictive truncation. With \code{zero.action = 0} (default), the function
#' automatically forbids bandwidths that yield erroneous values at the leave-one-out data point locations only. With \code{zero.action = 1}, the minimum machine value (see \code{.Machine$double.xmin} at the prompt) is
#' used to replace these individual leave-one-out values. When \code{zero.action = 2}, the minimum value of the valid (greater than zero) leave-one-out values is used to replace any erroneous leave-one-out values.
#' 
#' 
#' 
#' @aliases LIK.density
#' 
#' @rdname CV
#' 
#' @param pp An object of class \code{\link[spatstat.geom]{ppp}} giving the observed
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
#'   the number of CPU cores to use when \code{!auto.optim} \bold{and} \code{type = "fixed"}. Experimental. Test
#'   your system first using \code{parallel::detectCores()} to identify the
#'   number of cores available to you.
#' @param verbose Logical value indicating whether to provide function progress
#'   commentary.
#' @param type A character string; \code{"fixed"} (default) performs classical leave-one-out
#'   cross-validation for the fixed-bandwidth estimator. Alternatively, \code{"adaptive"} utilises
#'   multiscale adaptive kernel estimation (Davies & Baddeley, 2018) to run the cross-validation
#'   in an effort to find a suitable global bandwidth for the adaptive estimator. Note that data points are not `left out' of
#'   the pilot density estimate when using this option (this capability is currently in development). See also the entry for \code{...}.
#' @param zero.action A numeric integer, either \code{-1}, \code{0} (default), \code{1} or \code{2} controlling how the function should behave in response to numerical errors at very small bandwidths, when such a bandwidth results in one or more zero or negative density values during the leave-one-out computations. See `Details'.
#' @param ... Additional arguments controlling pilot density estimation and multi-scale bandwidth-axis
#'   resolution when \code{type = "adaptive"}. Relevant arguments are \code{hp}, \code{pilot.density},
#'   \code{gamma.scale}, and \code{trim} from \code{\link{bivariate.density}}; and \code{dimz} from 
#'   \code{\link{multiscale.density}}. If \code{hp} is missing and required, the function makes a (possibly recursive)
#'   call to \code{LSCV.density} to set this using fixed-bandwidth LSCV. The remaining defaults are \code{pilot.density = pp},
#'   \code{gamma.scale = "geometric"}, \code{trim = 5}, and \code{dimz = resolution}.
#'   
#' @return A single numeric value of the estimated bandwidth (if
#'   \code{auto.optim = TRUE}). Otherwise, a \eqn{[}\code{seqres} \eqn{x} 2\eqn{]} matrix 
#'   giving the bandwidth sequence and corresponding CV
#'   function value.
#'
#' @section Warning: Leave-one-out CV for bandwidth selection in kernel
#' density estimation is notoriously unstable in practice and has a tendency to
#' produce rather small bandwidths, particularly for spatial data. Satisfactory bandwidths are not guaranteed
#' for every application; \code{zero.action} can curb adverse numeric effects for very small bandwidths during the optimisation procedures. This method can also be computationally expensive for
#' large data sets and fine evaluation grid resolutions. The user may also need to
#' experiment with adjusting \code{hlim} to find a suitable minimum.
#'
#' @author T. M. Davies
#'
#' @seealso \code{\link{SLIK.adapt}} and functions for bandwidth selection in package
#'   \code{\link{spatstat}}: \code{\link[spatstat.core]{bw.diggle}};
#'   \code{\link[spatstat.core]{bw.ppl}}; \code{\link[spatstat.core]{bw.scott}};
#'   \code{\link[spatstat.core]{bw.frac}}.
#'
#' @references
#' 
#' Davies, T.M. and Baddeley A. (2018), Fast computation of
#' spatially adaptive kernel estimates, \emph{Statistics and Computing}, \bold{28}(4), 937-956.
#' 
#' Silverman, B.W. (1986), \emph{Density Estimation for Statistics
#' and Data Analysis}, Chapman & Hall, New York.
#'
#' Wand, M.P. and Jones,
#' C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.
#'
#' @examples
#' 
#' data(pbc)
#' pbccas <- split(pbc)$case
#' 
#' LIK.density(pbccas)
#' LSCV.density(pbccas)
#'
#' \donttest{
#' #* FIXED 
#' 
#' # custom limits
#' LIK.density(pbccas,hlim=c(0.01,4))
#' LSCV.density(pbccas,hlim=c(0.01,4))
#' 
#' # disable edge correction
#' LIK.density(pbccas,hlim=c(0.01,4),edge=FALSE)
#' LSCV.density(pbccas,hlim=c(0.01,4),edge=FALSE)
#' 
#' # obtain objective function
#' hcv <- LIK.density(pbccas,hlim=c(0.01,4),auto.optim=FALSE)
#' plot(hcv);abline(v=hcv[which.max(hcv[,2]),1],lty=2,col=2)
#' 
#' #* ADAPTIVE
#' LIK.density(pbccas,type="adaptive")
#' LSCV.density(pbccas,type="adaptive")
#'  
#' # change pilot bandwidth used
#' LIK.density(pbccas,type="adaptive",hp=2)
#' LSCV.density(pbccas,type="adaptive",hp=2)
#' } 
#' 
#' @export
LSCV.density <- function(pp,hlim=NULL,hseq=NULL,resolution=64,edge=TRUE,auto.optim=TRUE,
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
      result <- suppressWarnings(optimise(LSCV.density.spatial.single,interval=hlim,pp=pp,res=resolution,edge=edge,za=zero.action)$minimum)
      if(verbose) cat("Done.\n")
    } else {
      if(is.null(hseq)) hseq <- seq(hlim[1],hlim[2],length=seqres)
      hn <- length(hseq)
      if(is.null(parallelise)){
        lscv.vec <- rep(NA,hn)
        if(verbose) pb <- txtProgressBar(1,hn)
        for(i in 1:hn){
          lscv.vec[i] <- LSCV.density.spatial.single(hseq[i],pp,resolution,edge,za=zero.action)
          if(verbose) setTxtProgressBar(pb,i)
        }
        if(verbose) close(pb)
      } else {
        ncores <- detectCores()
        if(verbose) cat(paste("Evaluating criterion on",parallelise,"/",ncores,"cores..."))
        if(parallelise>ncores) stop("cores requested exceeds available count")
        registerDoParallel(cores=parallelise)
        lscv.vec <- foreach(i=1:hn,.packages="spatstat",.combine=c) %dopar% {
          return(LSCV.density.spatial.single(hseq[i],pp,resolution,edge,zero.action))
        }
        if(verbose) cat("Done.\n")
      }
      result <- cbind(hseq,lscv.vec)
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
    msobject <- multiscale.density(pp,h0=hhash,hp=hp,h0fac=hlim/hhash,edge=ifelse(edge,"uniform","none"),resolution=resolution,dimz=dimz,gamma.scale=gamma.scale,trim=trim,intensity=TRUE,pilot.density=pilot.density,verbose=FALSE)
    if(verbose) cat("Done.\n")
    
    h0range <- range(as.numeric(names(msobject$z)))
    if(auto.optim){
      if(verbose) cat("Searching for optimal h0 in ",prange(h0range),"...",sep="")
      h0opt <- suppressWarnings(optimise(ms.loo,interval=h0range,object=msobject,za=zero.action)$minimum)
      if(verbose) cat("Done.\n")
      return(h0opt)
    } else {
      if(is.null(hseq)) hseq <- seq(h0range[1],h0range[2],length=seqres)
      hn <- length(hseq)
      lscv.vec <- rep(NA,hn)
      if(verbose) pb <- txtProgressBar(1,hn)
      for(i in 1:hn){
        lscv.vec[i] <- ms.loo(hseq[i],msobject,za=zero.action)
        if(verbose) setTxtProgressBar(pb,i)
      }
      if(verbose) close(pb)
      
      result <- cbind(hseq,lscv.vec)
      dimnames(result)[[2]] <- c("h0","CV")
    }
  } else stop("invalid 'type'")
  
  return(result)
}
