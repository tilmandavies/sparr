#' Jointly optimal fixed bandwidth for the spatial relative risk function
#' 
#' Jointly optimal, common case-control isotropic fixed bandwidth for use in
#' the kernel-smoothed relative risk function.
#' 
#' Given the established preference of using a common bandwidth for both case
#' and control density estimates when constructing a fixed bandwidth relative
#' risk surface, This function calculates a `jointly optimal', common isotropic
#' LSCV bandwidth for the (Gaussian) kernel-smoothed relative risk function
#' (case-control density-ratio). It can be shown that choosing a bandwidth that
#' is equal for both case and control density estimates is preferable to
#' computing `separately optimal' bandwidths (Kelsall and Diggle, 1995).
#' \itemize{
#'   \item\code{method = "kelsall-diggle"}: the function computes the
#'     common bandwidth which minimises the approximate mean integrated squared
#'     error (MISE) of the log-transformed risk surface (Kelsall and Diggle, 1995).
#'   \item\code{method = "hazelton"}: the function minimises a
#'     \emph{weighted-by-control} MISE of the (raw) relative risk function
#'     (Hazelton, 2008).
#'   \item\code{method = "davies"}: the optimal bandwidth is
#'     one that minimises a crude plug-in approximation to the \emph{asymptotic}
#'     MISE (Davies, 2013).
#' }
#' 
#' @param f Either a pre-calculated object of class \code{\link{bivden}}
#'   representing the `case' (numerator) density estimate, or an object of class
#'   \code{\link[spatstat]{ppp}} giving the observed case data. Alternatively, if
#'   \code{f} is \code{\link[spatstat]{ppp}} object with dichotomous
#'   factor-valued \code{\link[spatstat]{marks}}, the function treats the first
#'   level as the case data, and the second as the control data, obviating the
#'   need to supply \code{g}.
#' @param g As for \code{f}, for the `control' (denominator) density; this
#'   object must be of the same class as \code{f}. Ignored if, as stated above,
#'   \code{f} contains both case and control observations.
#' @param hlim An optional vector of length 2 giving the limits of the
#'   optimisation routine with respect to the bandwidth. If unspecified, the
#'   function attempts to choose this automatically.
#' @param hseq An optional increasing sequence of bandwidth values at which to
#'   manually evaluate the optimisation criterion. Used only in the case
#'   \code{(!auto.optim && is.null(hlim))}.
#' @param method A character string controlling the selector to use. There are
#'   three types, based on either the mean integrated squared error (MISE)
#'   (Kelsall and Diggle, 1995; default -- \code{method = "kelsall-diggle"}); a
#'   weighted MISE (Hazelton, 2008 -- \code{method = "hazelton"}); or an
#'   approximation to the asymptotic MISE (Davies, 2013 -- \code{method =
#'   "davies"}). See `Details'.
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
#'
#' @return A single numeric value of the estimated bandwidth (if
#'   \code{auto.optim = TRUE}). Otherwise, a list of two numeric vectors of equal
#'   length giving the bandwidth sequence (as \code{hs}) and corresponding CV
#'   function value (as \code{CV}).
#'
#' @section Warning: The jointly optimal bandwidth selector can be
#' computationally expensive for large data sets and fine evaluation grid
#' resolutions. The user may need to experiment with adjusting \code{hlim} to
#' find a suitable minimum.
#'
#' @author T. M. Davies
#'
#' @seealso \code{\link{bivariate.density}}
#'
#' @references
#' Davies, T. M. (2013), Jointly optimal bandwidth selection for
#' the planar kernel-smoothed density-ratio, \emph{Spatial and Spatio-temporal
#' Epidemiology}, \bold{5}, 51-65.
#'
#' Hazelton, M. L. (2008), Letter to the
#' editor: Kernel estimation of risk surfaces without the need for edge
#' correction, \emph{Statistics in Medicine}, \bold{27}, 2269-2272.
#'
#' Kelsall, J.E. and Diggle, P.J. (1995), Kernel estimation of relative risk,
#' \emph{Bernoulli}, \bold{1}, 3-16.
#'
#' Silverman, B.W. (1986), \emph{Density Estimation for Statistics and Data Analysis},
#' Chapman & Hall, New York.
#'
#' Wand, M.P. and Jones, C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall,
#' London.
#'
#' @examples
#' 
#' ## To be filled
#' 
#' @export
LSCV.risk <- function(f, g = NULL, hlim = NULL, hseq = NULL, method = c("kelsall-diggle", "hazelton", "davies"),
                      resolution = 64, edge = TRUE, auto.optim = TRUE, seqres = 30,
                      parallelise = NULL, verbose = TRUE){
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
      marks(f) <- NULL
      marks(g) <- NULL
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
      marks(f) <- NULL
      marks(g) <- NULL
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
