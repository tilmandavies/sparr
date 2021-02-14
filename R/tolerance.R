#' Tolerance by \emph{p}-value surfaces
#' 
#' Calculates a \emph{p}-value surface based on asymptotic theory or
#' Monte-Carlo (MC) permutations describing the extremity of risk given a fixed
#' or adaptive kernel-smoothed density-ratio, allowing the drawing of
#' \emph{tolerance contours}.
#' 
#' This function implements developments in Hazelton and Davies (2009) (fixed)
#' and Davies and Hazelton (2010) (adaptive) to compute pointwise
#' \emph{p}-value surfaces based on asymptotic theory of kernel-smoothed
#' relative risk surfaces. Alternatively, the user may elect to calculate the
#' \emph{p}-value surfaces using Monte-Carlo methods (see Kelsall and Diggle,
#' 1995). Superimposition upon a plot of the risk surface contours of these
#' \emph{p}-values at given significance levels (i.e. ``tolerance contours'')
#' can be an informative way of exploring the statistical significance of the
#' extremity of risk across the defined study region.
#' 
#' Implementation of the Monte-Carlo method simply involves random allocation of case/control marks and
#' re-estimation of the risk surface \code{ITER} times, against which the
#' original estimate is compared.  While not dependent on asymptotic theory, it is
#' computationally expensive, and it has been suggested that it might have some
#' undesirable practical consequences in certain settings (Hazelton and Davies,
#' 2009). When performing the MC simulations, the same global (and pilot, if
#' necessary) bandwidths and edge-correction regimens are employed as were used
#' in the initial density estimates of the observed data. With regard to
#' arguments to be passed to internal calls of \code{\link{risk}}, the user
#' should take care to use \code{...} to set the \code{epsilon} value to match
#' that which was used in creation of the object passed to \code{rs} (if this
#' was set to a non-default value). Furthermore, if performing MC simulations
#' for the adaptive relative risk function, the function borrows the value of
#' the \code{beta} argument to speed things up via partitioning, and the user
#' should additionally access \code{...} to set the same \code{pilot.symmetry}
#' value as was used for creation of the object passed to \code{rs}, in the
#' same way as for any non-default use of \code{epsilon}. This will ensure the
#' simulations are all performed under the same conditions as were used to estimate the original risk
#' function.
#' 
#' 
#' 
#' 
#' @param rs An object of class \code{\link{rrs}} giving the estimated relative
#'   risk function for which to calculate the \emph{p}-value surface.
#' @param method A character string specifying the method of calculation.
#'   \code{"ASY"} (default) instructs the function to compute the \emph{p}-values
#'   using asymptotic theory. \code{"MC"} computes the values by random
#'   permutations of the data. See `Details'.
#' @param ref.density Required if \code{rs} is based on fixed-bandwidth
#'   estimates of the case and control densities and \code{method = "ASY"}.
#'   Either a pixel \code{\link[spatstat.geom]{im}}age or an object of class
#'   \code{\link{bivden}} giving the reference density to use in asymptotic
#'   formulae. May be unnormalised. Ignored if \code{rs} is based on adaptive
#'   kernel estimates or if \code{method = "MC"}.
#' @param beta A numeric value \eqn{0 <} \code{beta} \eqn{< 1} giving the
#'   fineness of the adaptive bandwidth partitioning to use for calculation of
#'   the required quantities for asymptotic adaptive \emph{p}-value surfaces.
#'   Smaller values provide more accurate bandwidth bins at the cost of
#'   additional computing time, see Davies and Baddeley (2018); the default is
#'   sensible in most cases. Ignored if \code{rs} is based on fixed-bandwidth
#'   kernel estimates.
#' @param ITER Number of iterations for the Monte-Carlo permutations. Ignored
#'   if \code{method = "ASY"}.
#' @param parallelise Numeric argument to invoke parallel processing, giving
#'   the number of CPU cores to use when \code{method = "MC"}. Experimental. Test
#'   your system first using \code{parallel::detectCores()} to identify the
#'   number of cores available to you.
#' @param verbose Logical value indicating whether to print function progress
#'   during execution.
#' @param ...  Additional arguments to be passed to \code{\link{risk}} when
#'   \code{method = "MC"}. While most information needed for the MC repetitions
#'   is implicitly gleaned from the object passed to \code{rs}, this ellipsis is
#'   typically used to set the appropriate \code{epsilon} and
#'   \code{pilot.symmetry} values for the internal calls to \code{\link{risk}}.
#'
#' @return A pixel \code{\link[spatstat.geom]{im}}age of the estimated
#' \emph{p}-value surface.
#'
#' @note The returned \emph{p}-values are geared so that ``smallness''
#' corresponds to statistical significance of elevated risk, that is, an
#' upper-tailed test. The complement of the \emph{p}-values will yeild
#' significance of reduced risk; a lower-tailed test. When using
#' \code{\link{tol.contour}}, the user can control what type of contours to
#' display.
#'
#' @author T. M. Davies
#'
#' @references
#'
#' Davies, T.M. and Baddeley A. (2018), Fast computation of
#' spatially adaptive kernel estimates, \emph{Statistics and Computing}, \bold{28}(4), 937-956.
#'
#' Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative
#' risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.
#' 
#' Davies, T.M., Jones, K. and Hazelton, M.L. (2016), Symmetric adaptive smoothing regimens for estimation of the spatial
#' relative risk function, \emph{Computational Statistics & Data Analysis},
#' \bold{101}, 12-28.
#'
#' Hazelton, M.L. and Davies, T.M. (2009), Inference based on kernel estimates
#' of the relative risk function in geographical epidemiology,
#' \emph{Biometrical Journal}, \bold{51}(1), 98-109.
#'
#' Kelsall, J.E. and Diggle, P.J. (1995), Kernel estimation of relative risk, \emph{Bernoulli},
#' \bold{1}, 3-16.
#'
#' @examples
#' 
#' \donttest{
#' 
#' data(pbc)
#' h0 <- LSCV.risk(pbc,method="hazelton");h0
#' pbccas <- split(pbc)[[1]]
#' pbccon <- split(pbc)[[2]]
#' 
#' # ASY
#' riskfix <- risk(pbc,h0=h0)
#' fixtol1 <- tolerance(riskfix,ref.density=density(pbc,OS(pbc)))
#' 
#' riskada <- risk(pbc,h0=h0,adapt=TRUE,hp=NS(pbc),pilot.symmetry="pooled",davies.baddeley=0.025)
#' adatol1 <- tolerance(riskada)
#' 
#' par(mfrow=c(1,2))
#' plot(riskfix)
#' tol.contour(fixtol1,levels=c(0.1,0.05,0.01),lty=3:1,add=TRUE)
#' plot(riskada)
#' tol.contour(adatol1,levels=c(0.1,0.05,0.01),lty=3:1,add=TRUE)
#' 
#'
#' # MC
#' fixtol2 <- tolerance(riskfix,method="MC",ITER=200) 
#' adatol2 <- tolerance(riskada,method="MC",ITER=200,parallelise=4) # ~1 minute with parallelisation
#' par(mfrow=c(1,2))
#' plot(riskfix)
#' tol.contour(fixtol2,levels=c(0.1,0.05,0.01),lty=3:1,add=TRUE)
#' plot(riskada)
#' tol.contour(adatol2,levels=c(0.1,0.05,0.01),lty=3:1,add=TRUE)
#' }
#' 
#' 
#' @export
tolerance <- function(rs, method = c("ASY", "MC"), ref.density = NULL, beta = 0.025,
                      ITER = 100, parallelise = NULL, verbose = TRUE, ...){
  if(!inherits(rs,"rrs")) stop("'rs' argument must be of class \"rrs\"")
  
  meth <- method[1]
  adaf <- !is.null(rs$f$him)
  adag <- !is.null(rs$g$him)
  ada <- adaf + adag
  if(ada==1) stop("case/control smoothing regimens (fixed or adaptive) must be identical")
  
  if(meth=="ASY"){
    if(ada==2){
      psurf <- tol.asy.ada(rs$f,rs$g,beta,verbose)$p
    } else {
      if(is.null(ref.density)) stop("must supply 'ref.density' for fixed-bandwidth asymptotic tolerance contours")
      if((!inherits(ref.density,"bivden"))&&(!inherits(ref.density,"im"))) stop("'ref.density' must be of class \"bivden\" or \"im\"")
      if(is.im(ref.density)){
        ref.density <- list(z=ref.density,q=NULL)
      #  was:
      #   rq <- ref.density
      # 	rq$v[!is.na(rq$v)] <- 1
      # 	ref.density <- list(z=ref.density,q=rq)
      }
      ref.density$z <- ref.density$z/integral(ref.density$z)
      if(!compatible(rs$f$z,rs$g$z,ref.density$z)) stop("incompatible 'ref.density'... must be evaluated on domain identical to case/control densities")
      psurf <- tol.asy.fix(rs$f,rs$g,ref.density,verbose)$p
    }
  } else if(meth=="MC"){
    ITER <- checkit(ITER,"'ITER'")
    if(ada==2){
      psurf <- im(tol.mc.ada(rs,round(ITER),parallelise,verbose,davies.baddeley=beta,...),xcol=rs$rr$xcol,yrow=rs$rr$yrow)
    } else {
      psurf <- im(tol.mc.fix(rs,round(ITER),parallelise,verbose,...),xcol=rs$rr$xcol,yrow=rs$rr$yrow)
    }
  } else stop("invalid 'method' argument")
  
  return(psurf)
}
  
