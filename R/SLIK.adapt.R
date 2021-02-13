#' Simultaneous global/pilot likelihood bandwidth selection
#' 
#' Isotropic global and pilot bandwidth selection for adaptive density/intensity
#' based on likelihood cross-validation.
#' 
#' This function is a generalisation of \code{\link{LIK.density}}, and is used in attempts to simultaneously choose
#' an optimal global and pilot bandwidth for adaptive kernel density estimates. Where \code{\link{LIK.density}} for adaptive
#' estimates assumes the pilot density is held constant (and is not subject to the leave-one-out operations), this function
#' allows the pilot bandwidth to vary alongside the global.
#' 
#' Thus, in contrast to \code{\link{LIK.density}} the internal leave-one-out operations now also affect the
#' pilot estimation stage. Hence, the set of variable bandwidths changes as each point is left out. In turn, this means the leave-one-out operations must
#' be computed by brute force, and this is computationally expensive.
#' 
#' Identifiability problems can sometimes arise when the global and pilot bandwidths are allowed to `float freely' in the bivariate optimisation routine, which is the default
#' behaviour of the function (with \code{hold = FALSE}). This can be curbed by setting \code{hold = TRUE}, which forces both the global and pilot
#' to be held at the same value during optimisation. Doing this also has the beneficial side effect of turning the problem into one of univariate optimisation, thereby reducing total computational cost. Current work (Davies & Lawson, 2018) provides some empirical evidence that this strategy performs quite well in practice.
#' 
#' Like \code{\link{LSCV.density}} and \code{\link{LIK.density}}, the argument \code{zero.action} can be used to control the level of severity in response to small bandwidths that result (due to numerical error) in at least one density value being zero or less. 
#' When this argument is passed a vector of length 2, the first entry corresponds to the global bandwidth (and hence refers to checks of the final adaptive density estimate and its leave-one-out values) and the second to the pilot bandwidth (and hence checks the fixed-bandwidth pilot density and its leave-one-out values).
#' Alternatively a single value may be supplied, which will be taken to be the same for both global and pilot.
#' See the help page for \code{\link{LIK.density}} for an explanation of the four allowable values (\code{-1}, \code{0}, \code{1}, \code{2}) for each component of this argument.
#' 
#' 
#' 
#' @rdname SLIK.adapt
#' 
#' @param pp An object of class \code{\link[spatstat.geom]{ppp}} giving the observed
#'   2D data to be smoothed.
#' @param hold Logical value indicating whether to hold the global and pilot bandwidths equal throughout the
#'   optimisation; defaults to \code{TRUE}. See `Details'. 
#' @param hlim An optional vector of length 2 giving the limits of the
#'   optimisation routine with respect to the bandwidth when \code{hold = TRUE}. If unspecified, the
#'   function attempts to choose this automatically. Ignored when \code{hold = FALSE}.
#' @param start A positively-valued numeric vector of length 2 giving the starting values to be used for the global/pilot
#'   optimisation routine when \code{hold = FALSE}. Defaults to the oversmoothing bandwidth (\code{\link{OS}}) for both values;
#'   ignored when \code{hold = TRUE}.
#' @param edge Logical value indicating whether to edge-correct the density
#'   estimates used.
#' @param parallelise Numeric argument to invoke parallel processing in the brute force leave-one-out calculations, giving
#'   the number of CPU cores to use. Experimental. Test
#'   your system first using \code{parallel::detectCores()} to identify the
#'   number of cores available to you. If \code{NA} (default), no parallelisation performed and a single loop is used.
#' @param verbose Logical value indicating whether to provide function progress
#'   commentary.
#' @param zero.action A numeric vector of length 2, each value being either \code{-1}, \code{0} (default), \code{1} or \code{2} controlling how the function should behave in response to numerical errors at very small bandwidths, when such a bandwidth results in one or more zero or negative density values during the leave-one-out computations. See `Details'.
#' @param optim.control An optional list to be passed to the \code{control} argument of \code{\link[stats]{optim}} for further control over the numeric optimisation when \code{hold = FALSE}. See the documentation for \code{\link[stats]{optim}} for further details.
#' @param ... Additional arguments controlling density estimation for the internal calculations. Relevant arguments are \code{resolution}, \code{gamma.scale}, and \code{trim}. If unsupplied these default to \code{64}, \code{"geometric"}, and \code{5} respectively; see \code{\link{bivariate.density}} for a further explanation of these arguments.
#'   
#' @return A numeric vector of length 2 giving the likelihood-maximised global and pilot bandwidths.
#' 
#' @section Note: While theoretically valid, this is a largely experimental function. There is presently little in the literature to suggest how well this
#' type of simultaneous global/pilot bandwidth selection might perform in practice. Current research efforts (Davies & Lawson, 2018)
#' seek in part to address these questions.
#' 
#' 
#'
#' @author T. M. Davies
#'
#' @seealso Functions for bandwidth selection in package
#'   \code{\link{spatstat}}: \code{\link[spatstat.core]{bw.diggle}};
#'   \code{\link[spatstat.core]{bw.ppl}}; \code{\link[spatstat.core]{bw.scott}};
#'   \code{\link[spatstat.core]{bw.frac}}.
#'
#' @references
#' Davies, T.M. and Lawson, A.B. (2018), An evaluation of likelihood-based bandwidth selectors for spatial and spatiotemporal kernel estimates, \emph{Submitted for publication}.
#' 
#' Silverman, B.W. (1986), \emph{Density Estimation for Statistics
#' and Data Analysis}, Chapman & Hall, New York.
#'
#' Wand, M.P. and Jones,
#' C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.
#'
#' @examples
#' 
#'
#' \donttest{
#' 
#' data(pbc)
#' pbccas <- split(pbc)$case
#' 
#' SLIK.adapt(pbccas)
#' SLIK.adapt(pbccas,hold=TRUE)
#' 
#' } 
#' 
#' @export
SLIK.adapt <- function(pp,hold=TRUE,start=rep(OS(pp),2),hlim=NULL,edge=TRUE,zero.action=c(-1,0),optim.control=list(),parallelise=NULL,verbose=TRUE,...){
  if(class(pp)!="ppp") stop("data object 'pp' must be of class \"ppp\"")
  W <- Window(pp)
  # resolution <- checkit(resolution,"'resolution'")
  
  # Handle the dots #
  ellip <- list(...)
  
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
  
  if(is.null(ellip$resolution)){
    resolution <- 64
  } else {
    resolution <- ellip$resolution
  }
  ##
  
  if(length(zero.action)==1) zero.action <- rep(zero.action,2)
  zero.action <- zero.action[1:2]
  if((!zero.action[1]%in%c((-1):2))||(!zero.action[2]%in%c((-1):2))) stop("invalid 'zero.action'")
  
  # optimise/optim #
  if(!hold){
    result <- optim(start,loowrap_nohold,pp=pp,edge=edge,gamma.scale=gamma.scale,trim=trim,resolution=resolution,vbs=verbose,parallel=parallelise,za=zero.action,control=optim.control)$par
    result <- as.numeric(result)
  } else {
    if(is.null(hlim)){
      ppu <- pp
      marks(ppu) <- NULL
      md <- min(nndist(unique(ppu)))
      hlim <- c(md,max(md*50,min(diff(W$xrange),diff(W$yrange))/6))
    } else {
      hlim <- checkran(hlim,"'hlim'")
    }
    result <- suppressWarnings(optimise(loowrap_hold,interval=hlim,pp=pp,edge=edge,gamma.scale=gamma.scale,trim=trim,resolution=resolution,vbs=verbose,parallel=parallelise,za=zero.action)$minimum)
    result <- as.numeric(rep(result,2))
  }
  names(result) <- c("h0","hp")
  return(result)
}

loowrap_hold <- function(hh,...) return(loowrap_nohold(c(hh,hh),...))
  
loowrap_nohold <- function(h0hp,pp,edge,gamma.scale,trim,resolution,vbs,parallel,za){
  W <- Window(pp)
  
  if(any(h0hp<=0)||any(h0hp>100*max(c(diff(W$xrange),diff(W$yrange))))) return(NA)
  
  if(vbs) cat("h0: ",h0hp[1],"; hp: ",h0hp[2],"\n",sep="")
  
  if(za[2]==-1){
    stopper <- density(pp,sigma=h0hp[2],edge=edge,dimyx=64,weights=NULL)
    if(any(stopper<=0)) return(Inf)
  }
  if(za[1]==-1){
    stopper <- bivariate.density(pp,h0hp[1],h0hp[2],adapt=TRUE,edge=ifelse(edge,"uniform","none"),trim=trim,gamma.scale=gamma.scale,verbose=FALSE,resolution=64,weights=NULL,davies.baddeley=0.05)$z
    if(any(stopper<=0)) return(Inf)
  }

  loovals <- bivden.LOO(pp,h0hp[1],h0hp[2],edge=edge,trim=trim,gamma.scale=gamma.scale,resolution=resolution,parallel=parallel,weights=NULL,za[2])[[1]]
  
  if(any(loovals<=0)){
    if(za[1]==0){
      loovals[loovals<=0] <- NA
    } else if(za[1]==1){
      loovals <- posifybivden(loovals)
    } else if(za[1]==2){
      loovals[loovals<=0] <- min(loovals[loovals>0])
    }
  }
  
  return(-mean(log(loovals)))
}



  
