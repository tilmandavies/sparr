#' Jointly optimal bandwidth selection for the spatial relative risk function
#' 
#' Methods to find a jointly optimal, common case-control isotropic bandwidth for use in
#' estimation of the fixed or adaptive kernel-smoothed relative risk function.
#' 
#' Given the established preference of using a common bandwidth for both case
#' and control density estimates when constructing a relative
#' risk surface, This function calculates a `jointly optimal', common isotropic
#' LSCV bandwidth for the (Gaussian) kernel-smoothed relative risk function
#' (case-control density-ratio). It can be shown that choosing a bandwidth that
#' is equal for both case and control density estimates is preferable to
#' computing `separately optimal' bandwidths (Kelsall and Diggle, 1995). The user
#' can choose to either calculate a common smoothing parameter for a fixed-bandwidth
#' relative risk surface (\code{type = "fixed"}; default), or a common global bandwidth for
#' an adaptive risk surface (\code{type = "adaptive"}). See further comments below.
#' 
#' 
#' 
#' \itemize{
#'   \item\code{method = "kelsall-diggle"}: the function computes the
#'     common bandwidth which minimises the approximate mean integrated squared
#'     error (MISE) of the log-transformed risk surface (Kelsall and Diggle, 1995).
#'   \item\code{method = "hazelton"}: the function minimises a
#'     \emph{weighted-by-control} MISE of the (raw) relative risk function
#'     (Hazelton, 2008).
#'   \item\code{method = "davies"}: the optimal bandwidth is
#'     one that minimises a crude plug-in approximation to the \emph{asymptotic}
#'     MISE (Davies, 2013). Only possible for \code{type = "fixed"}.
#' }
#' 
#' For jointly optimal, common global bandwidth selection when \code{type = "adaptive"}, the
#' optimisation routine utilises \code{\link{multiscale.density}}. Like \code{\link{LSCV.density}},
#' the leave-one-out procedure does not affect the pilot density, for which additional
#' control is offered via the \code{hp} and \code{pilot.symmetry} arguments. The user has the option of
#' obtaining a so-called \emph{symmetric} estimate (Davies et al. 2016) via
#' \code{pilot.symmetry}. This amounts to choosing the same pilot density for
#' both case and control densities. By choosing \code{"none"} (default), the
#' result uses the case and control data separately for the fixed-bandwidth
#' pilots, providing the original asymmetric density-ratio of Davies and
#' Hazelton (2010). By selecting either of \code{"f"}, \code{"g"}, or
#' \code{"pooled"}, the pilot density is calculated based on the case, control,
#' or pooled case/control data respectively (using \code{hp[1]} as the fixed
#' bandwidth). Davies et al. (2016) noted some beneficial practical behaviour
#' of the symmetric adaptive surface over the asymmetric. (The pilot bandwidth(s), if not supplied in \code{hp}, are calculated
#' internally via default use of \code{\link{LSCV.density}}, using the requested symmetric-based data set, or separately with respect to the case and control datasets \code{f} and \code{g} if
#' \code{pilot.symmetry = "none"}.)
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
#' @param type A character string; \code{"fixed"} (default) performs classical leave-one-out
#'   cross-validation for a jointly optimal fixed bandwidth. Alternatively, \code{"adaptive"} utilises
#'   multiscale adaptive kernel estimation (Davies & Baddeley, 2018) to run the cross-validation
#'   in an effort to find a suitable jointly optimal, common global bandwidth for the adaptive relative risk function. See `Details'.
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
#' @param hp A single numeric value or a vector of length 2 giving the pilot
#'   bandwidth(s) to be used for estimation of the pilot
#'   densities for adaptive risk surfaces. Ignored if \code{type = "fixed"}.
#' @param pilot.symmetry A character string used to control the type of
#'   symmetry, if any, to use for the bandwidth factors when computing an
#'   adaptive relative risk surface. See `Details'. Ignored if \code{type = "fixed"}.
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
#' @param ... Additional arguments such as \code{dimz} and \code{trim} to be passed to
#'   the internal calls to \code{\link{multiscale.density}}.
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
#' 
#' Davies, T. M. (2013), Jointly optimal bandwidth selection for
#' the planar kernel-smoothed density-ratio, \emph{Spatial and Spatio-temporal
#' Epidemiology}, \bold{5}, 51-65.
#'
#' Davies, T.M. and Baddeley A. (2018), Fast computation of
#' spatially adaptive kernel estimates, \emph{Statistics and Computing}, \bold{28}(4), 937-956.
#' 
#' Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel
#' estimation of spatial relative risk, \emph{Statistics in Medicine},
#' \bold{29}(23) 2423-2437.
#' 
#' Davies, T.M., Jones, K. and Hazelton, M.L.
#' (2016), Symmetric adaptive smoothing regimens for estimation of the spatial
#' relative risk function, \emph{Computational Statistics & Data Analysis},
#' \bold{101}, 12-28.
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
#' \donttest{
#' 
#' data(pbc)
#' pbccas <- split(pbc)$case
#' pbccon <- split(pbc)$control
#' 
#' # FIXED (for common h)
#' 
#' LSCV.risk(pbccas,pbccon)
#' LSCV.risk(pbccas,pbccon,method="hazelton")
#' hcv <- LSCV.risk(pbccas,pbccon,method="davies",auto.optim=FALSE)
#' plot(hcv[,1],log(hcv[,2]));abline(v=hcv[which.min(hcv[,2]),1],col=2,lty=2)
#' 
#' 
#' # ADAPTIVE (for common h0)
#'
#' LSCV.risk(pbccas,pbccon,type="adaptive")
#' 
#' # change pilot bandwidths used
#' LSCV.risk(pbccas,pbccon,type="adaptive",hp=c(OS(pbccas)/2,OS(pbccon)/2))
#' 
#' # specify pooled-data symmetric relative risk estimator 
#' LSCV.risk(pbccas,pbccon,type="adaptive",hp=OS(pbc),pilot.symmetry="pooled")
#' 
#' # as above, for Hazelton selector
#' LSCV.risk(pbccas,pbccon,type="adaptive",method="hazelton")
#' LSCV.risk(pbccas,pbccon,type="adaptive",method="hazelton",hp=c(OS(pbccas)/2,OS(pbccon)/2))
#' LSCV.risk(pbccas,pbccon,type="adaptive",method="hazelton",hp=OS(pbc),pilot.symmetry="pooled")
#' }
#' 
#' @export
LSCV.risk <- function(f, g = NULL, hlim = NULL, hseq = NULL, type = c("fixed", "adaptive"),
                      method = c("kelsall-diggle", "hazelton", "davies"),
                      resolution = 64, edge = TRUE, hp = NULL, pilot.symmetry = c("none","f","g","pooled"),
                      auto.optim = TRUE, seqres = 30,
                      parallelise = NA, verbose = TRUE, ...){
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
  typ <- type[1]
  if(meth=="davies"&&typ=="adaptive") stop("method = \"davies\" not possible for type = \"adaptive\"")

  if(typ=="fixed"){
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
        if(is.na(parallelise)){
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
        if(is.na(parallelise)){
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
        if(is.na(parallelise)){
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
      
  } else if(typ=="adaptive"){
    
    pilot.symmetry <- pilot.symmetry[1]
    pdat <- list()
    if(pilot.symmetry=="none"){
      pdat[[1]] <- f
      pdat[[2]] <- g
    } else if(pilot.symmetry=="f"){
      pdat[[1]] <- pdat[[2]] <- f
    } else if(pilot.symmetry=="g"){
      pdat[[1]] <- pdat[[2]] <- g
    } else if(pilot.symmetry=="pooled"){
      marks(f) <- NULL
      marks(g) <- NULL
      pooled <- suppressWarnings(superimpose(f,g))
      pdat[[1]] <- pdat[[2]] <- pooled
    } else {
      stop("invalid 'pilot.symmetry' argument")
    }
    
    if(!is.null(hp)){
      if(length(hp)>1){
        fp <- hp[1]
        gp <- hp[2]
      } else {
        fp <- gp <- hp[1]
      }
    } else {
      if(verbose) cat("Selecting pilot bandwidth(s)...")
      if(pilot.symmetry=="none"){
        if(verbose) cat("\n --f--\n")
        fp <- LSCV.density(f,verbose=FALSE)
        if(verbose) cat(" --g--\n")
        gp <- LSCV.density(g,verbose=FALSE)
      } else {
        fp <- gp <- LSCV.density(pdat[[1]],verbose=FALSE)
      }
      if(verbose) cat(paste("Done.\n   [ Using hp(f) =",fp,"\b; hp(g) =",gp,"]\n"))
    }
    
    hhash <- mean(hlim)
    if(verbose) cat("Computing multi-scale estimates...\n --f--\n")
    fms <- multiscale.density(f,h0=hhash,hp=fp,h0fac=hlim/hhash,edge=ifelse(edge,"uniform","none"),resolution=resolution,intensity=FALSE,pilot.density=pdat[[1]],verbose=FALSE,...)
    if(verbose) cat(" --g--\n")
    gms <- multiscale.density(g,h0=hhash,hp=gp,h0fac=hlim/hhash,edge=ifelse(edge,"uniform","none"),resolution=resolution,intensity=FALSE,pilot.density=pdat[[2]],verbose=FALSE,...)
    if(verbose) cat("Done.\n")
    
    h0range <- fms$h0range
    
    if(meth=="kelsall-diggle"){
      if(auto.optim){
        if(verbose) cat("Searching for optimal h0 in ",prange(h0range),"...",sep="")
        h0opt <- optimise(ms.loo.risk,interval=h0range,fob=fms,gob=gms,hazey=FALSE)$minimum
        if(verbose) cat("Done.\n")
        return(h0opt)
      } else {
        if(is.null(hseq)) hseq <- seq(h0range[1],h0range[2],length=seqres)
        hn <- length(hseq)
        if(is.na(parallelise)){
          lscv.vec <- rep(NA,hn)
          if(verbose) pb <- txtProgressBar(1,hn)
          for(i in 1:hn){
            lscv.vec[i] <- suppressWarnings(ms.loo.risk(hseq[i],fob=fms,gob=gms,hazey=FALSE))
            if(verbose) setTxtProgressBar(pb,i)
          }
          if(verbose) close(pb)
        } else {
          ncores <- detectCores()
          if(verbose) cat(paste("Evaluating criterion on",parallelise,"/",ncores,"cores..."))
          if(parallelise>ncores) stop("cores requested exceeds available count")
          registerDoParallel(cores=parallelise)
          lscv.vec <- foreach(i=1:hn,.combine=c) %dopar% { # .packages="spatstat"
            return(suppressWarnings(ms.loo.risk(hseq[i],fob=fms,gob=gms,hazey=FALSE)))
          }
          if(verbose) cat("Done.\n")
        }
      }
    } else if(meth=="hazelton"){
      if(auto.optim){
        if(verbose) cat("Searching for optimal h0 in ",prange(h0range),"...",sep="")
        h0opt <- optimise(ms.loo.risk,interval=h0range,fob=fms,gob=gms,hazey=TRUE)$minimum
        if(verbose) cat("Done.\n")
        return(h0opt)
      } else {
        if(is.null(hseq)) hseq <- seq(h0range[1],h0range[2],length=seqres)
        hn <- length(hseq)
        if(is.na(parallelise)){
          lscv.vec <- rep(NA,hn)
          if(verbose) pb <- txtProgressBar(1,hn)
          for(i in 1:hn){
            lscv.vec[i] <- suppressWarnings(ms.loo.risk(hseq[i],fob=fms,gob=gms,hazey=TRUE))
            if(verbose) setTxtProgressBar(pb,i)
          }
          if(verbose) close(pb)
        } else {
          ncores <- detectCores()
          if(verbose) cat(paste("Evaluating criterion on",parallelise,"/",ncores,"cores..."))
          if(parallelise>ncores) stop("cores requested exceeds available count")
          registerDoParallel(cores=parallelise)
          lscv.vec <- foreach(i=1:hn,.combine=c) %dopar% { # .packages="spatstat"
            return(suppressWarnings(ms.loo.risk(hseq[i],fob=fms,gob=gms,hazey=TRUE)))
          }
          if(verbose) cat("Done.\n")
        }
      }
    } else {
      stop("invalid 'method'")
    }
    
    result <- cbind(hseq,lscv.vec)
    dimnames(result)[[2]] <- c("h","CV")
    return(result)
    
  } else {
    stop("invalid 'type'")
  }

}


# if(!is.null(pilot.args)&&!is.list(pilot.args)) stop("'pilot.args' must be a list")
# 
# if(!is.null(pilot.args$pilot.density)){
#   fp <- gp <- NULL
#   if(is.list(pilot.args$pilot.density)&&!is.im(pilot.args$pilot.density)){
#     if(length(pilot.args$pilot.density)>1){
#       fpilot <- pilot.args$pilot.density[[1]]
#       gpilot <- pilot.args$pilot.density[[2]]
#       if(!is.im(fpilot)||!is.im(gpilot)) stop("pilot.args$pilot.density must be a pixel image ('spatstat' class 'im') or a list of two pixel images")
#     }
#   } else {
#     if(!is.im(pilot.args$pilot.density)||!is.ppp(pilot.args$pilot.density)) stop("pilot.args$pilot.density must be of class 'im' or 'ppp' (or a list of two)")
#     fpilot <- gpilot <- pilot.args$pilot.density
#   }
# } else {
#   fpilot <- gpilot <- NULL
#   
# if(!is.null(pilot.args$hp)){
#   if(length(pilot.args$hp)>1){
#     fp <- pilot.args$hp[1]
#     gp <- pilot.args$hp[2]
#   } else {
#     fp <- gp <- pilot.args$hp[1]
#   }
# } else {
#   if(verbose) cat("Selecting pilot bandwidths...\n --f--\n")
#   fp <- LSCV.density(f,verbose=FALSE)
#   if(verbose) cat(" --g--\n")
#   gp <- LSCV.density(g,verbose=FALSE)
#   if(verbose) cat(paste("Done.\n   [ Found hp(f) =",fp,"\b; hp(g) =",gp,"]\n"))
# }
# }

# if(!is.null(pilot.args$dimz)){
#   dimz <- pilot.args$dimz[1]
# } else {
#   dimz <- resolution
# }
# 
# if(!is.null(pilot.args$trim)){
#   trim <- pilot.args$trim
# } else {
#   trim <- 5
# }

#  
# the \code{pilot.args} argument. This should be supplied as a named
# list, with optional components \code{hp}, \code{pilot.density}, \code{dimz}, and \code{trim}.
# See the documentation for these arguments in \code{\link{multiscale.density}}. By default, \code{trim = 5};
# \code{dimz = resolution}; \code{pilot.density = NULL}; and 
# Otherwise, the \code{pilot.density} component can be a single
# pixel \code{\link[spatstat]{im}}age (defined on the same domain as the data in \code{f} and \code{g}, and also matching \code{resolution}), posing as the common pilot density (i.e. if the selected global bandwidth 
# is intended for a symmetric adaptive relative risk surface, see Davies et al. 2016). Aternatively, the \code{pilot.density} component can be provided as a list 
# of two pixel images -- one for the case density, the other for the control (in that order).
# The \code{hp} component is only used if \code{pilot.args$pilot.density} is unsupplied, in which case it should be a vector of length one or two giving either a common pilot bandwidth
# or the case and control pilot bandwidths respectively. Either way, unless \code{pilot.args$pilot.density} is a single pixel \code{\link[spatstat]{im}}age as noted above, 
# the pilot densities are computed separately using the case (\code{f}) and control (\code{g}) data supplied to the function for an asymmetric adaptive relative risk surface (Davies & Hazelton, 2010).
 
