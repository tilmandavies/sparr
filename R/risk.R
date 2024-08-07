#' Spatial relative risk/density ratio
#' 
#' Estimates a \emph{relative risk} function based on the ratio of two 2D
#' kernel density estimates.
#' 
#' The relative risk function is defined here as the ratio of the `case'
#' density to the `control' (Bithell, 1990; 1991). Using kernel density
#' estimation to model these densities (Diggle, 1985), we obtain a workable
#' estimate thereof. This function defines the risk function \emph{r} in the
#' following fashion: \cr\cr \emph{r}\code{ = (fd + epsilon*max(gd))/(gd +
#' epsilon*max(gd))}, \cr\cr where \code{fd} and \code{gd} denote the case and
#' control density estimates respectively. Note the (optional) additive
#' constants defined by \code{epsilon} times the maximum of each of the
#' densities in the numerator and denominator respectively (see Bowman and
#' Azzalini, 1997). A more recent shrinkage estimator developed by Hazelton (2023)
#' is also implemented.
#' 
#' The log-risk function \emph{rho}, given by \emph{rho} = log[\emph{r}], is
#' argued to be preferable in practice as it imparts a sense of symmetry in the
#' way the case and control densities are treated (Kelsall and Diggle,
#' 1995a;b). The option of log-transforming the returned risk function is
#' therefore selected by default.
#' 
#' When computing adaptive relative risk functions, the user has the option of
#' obtaining a so-called \emph{symmetric} estimate (Davies et al. 2016) via
#' \code{pilot.symmetry}. This amounts to choosing the same pilot density for
#' both case and control densities. By choosing \code{"none"} (default), the
#' result uses the case and control data separately for the fixed-bandwidth
#' pilots, providing the original asymmetric density-ratio of Davies and
#' Hazelton (2010). By selecting either of \code{"f"}, \code{"g"}, or
#' \code{"pooled"}, the pilot density is calculated based on the case, control,
#' or pooled case/control data respectively (using \code{hp[1]} as the fixed
#' bandwidth). Davies et al. (2016) noted some beneficial practical behaviour
#' of the symmetric adaptive surface over the asymmetric.
#' 
#' If the user selects \code{tolerate = TRUE}, the function internally computes
#' asymptotic tolerance contours as per Hazelton and Davies (2009) and Davies
#' and Hazelton (2010). When \code{adapt = FALSE}, the reference density
#' estimate (argument \code{ref.density} in \code{\link{tolerance}}) is taken
#' to be the estimated control density. The returned pixel
#' \code{\link[spatstat.geom]{im}}age of \emph{p}-values (see `Value') is
#' interpreted as an upper-tailed test i.e. smaller \emph{p}-values represent
#' greater evidence in favour of significantly increased risk. For greater
#' control over calculation of tolerance contours, use \code{\link{tolerance}}.
#' 
#' @aliases risk rrs
#'
#' @param f Either a pre-calculated object of class \code{\link{bivden}}
#'   representing the `case' (numerator) density estimate, or an object of class
#'   \code{\link[spatstat.geom]{ppp}} giving the observed case data. Alternatively, if
#'   \code{f} is \code{\link[spatstat.geom]{ppp}} object with dichotomous
#'   factor-valued \code{\link[spatstat.geom]{marks}}, the function treats the first
#'   level as the case data, and the second as the control data, obviating the
#'   need to supply \code{g}.
#' @param g As for \code{f}, for the `control' (denominator) density; this
#'   object must be of the same class as \code{f}. Ignored if, as stated above,
#'   \code{f} contains both case and control observations.
#' @param log Logical value indicating whether to return the (natural)
#'   log-transformed relative risk function as recommended by Kelsall and Diggle
#'   (1995a). Defaults to \code{TRUE}, with the alternative being the raw density
#'   ratio.
#' @param h0 A single positive numeric value or a vector of length 2 giving the
#'   global bandwidth(s) to be used for case/control density estimates;
#'   defaulting to a common oversmoothing bandwidth computed via \code{\link{OS}}
#'   on the pooled data using \code{nstar = "geometric"} if unsupplied. Ignored if \code{f} and \code{g} are
#'   already \code{\link{bivden}} objects.
#' @param hp A single numeric value or a vector of length 2 giving the pilot
#'   bandwidth(s) to be used for fixed-bandwidth estimation of the pilot
#'   densities for adaptive risk surfaces. Ignored if \code{adapt = FALSE} or if
#'   \code{f} and \code{g} are already \code{\link{bivden}} objects.
#' @param adapt A logical value indicating whether to employ adaptive smoothing
#'   for internally estimating the densities. Ignored if \code{f} and \code{g}
#'   are already \code{\link{bivden}} objects.
#' @param shrink A logical value indicating whether to compute the shrinkage estimator
#'   of Hazelton (2023). This is only possible for \code{adapt=FALSE}.
#' @param shrink.args A named list of optional arguments controlling the shrinkage estimator.
#'   Possible entries are \code{rescale} (a logical value indicating whether to integrate
#'   to one with respect to the control distribution over the window); \code{type} 
#'   (a character string stipulating the shrinkage methodology to be used, either the
#'   default \code{"lasso"} or the alternative \code{"Bithell"}); and \code{lambda}
#'   (a non-negative numeric value determining the degree of shrinkage towards uniform
#'   relative risk---when set to its default \code{NA}, it is selected via cross-validation).
#' @param tolerate A logical value indicating whether to internally calculate a
#'   corresponding asymptotic p-value surface (for tolerance contours) for the
#'   estimated relative risk function. See `Details'.
#' @param doplot Logical. If \code{TRUE}, an image plot of the estimated
#'   relative risk function is produced using various visual presets. If
#'   additionally \code{tolerate} was \code{TRUE}, asymptotic tolerance contours
#'   are automatically added to the plot at a significance level of 0.05 for
#'   elevated risk (for more flexible options for calculating and plotting
#'   tolerance contours, see \code{\link{tolerance}} and
#'   \code{\link{tol.contour}}).
#' @param pilot.symmetry A character string used to control the type of
#'   symmetry, if any, to use for the bandwidth factors when computing an
#'   adaptive relative risk surface. See `Details'. Ignored if \code{adapt =
#'   FALSE}.
#' @param epsilon A single non-negative numeric value used for optional scaling
#'   to produce additive constant to each density in the raw ratio (see
#'   `Details'). A zero value requests no additive constant (default).
#' @param verbose Logical value indicating whether to print function progress
#'   during execution.
#' @param ...  Additional arguments passed to any internal calls of
#'   \code{\link{bivariate.density}} for estimation of the requisite densities.
#'   Ignored if \code{f} and \code{g} are already \code{\link{bivden}} objects.
#'
#' @return An object of class \code{"rrs"}. This is a named list with the
#' following components:
#' \item{rr}{A pixel \code{\link[spatstat.geom]{im}}age of the
#'   estimated risk surface.}
#' \item{f}{An object of class \code{\link{bivden}}
#'   used as the numerator or `case' density estimate.}
#' \item{g}{An object of
#'   class \code{\link{bivden}} used as the denominator or `control' density
#'   estimate.}
#' \item{P}{Only included if \code{tolerate = TRUE}. A pixel
#'   \code{\link[spatstat.geom]{im}}age of the \emph{p}-value surface for tolerance
#'   contours; \code{NULL} otherwise.}
#'
#' @author T.M. Davies
#'
#' @references
#' Bithell, J.F. (1990), An application of density estimation to
#' geographical epidemiology, \emph{Statistics in Medicine}, \bold{9},
#' 691-701.
#'
#' Bithell, J.F. (1991), Estimation of relative risk functions,
#' \emph{Statistics in Medicine}, \bold{10}, 1745-1751.
#'
#' Bowman, A.W. and Azzalini A. (1997), \emph{Applied Smoothing Techniques for Data Analysis:
#' The Kernel Approach with S-Plus Illustrations}, Oxford University Press
#' Inc., New York.
#'
#' Davies, T.M. and Hazelton, M.L. (2010), Adaptive
#' kernel estimation of spatial relative risk, \emph{Statistics in Medicine},
#' \bold{29}(23) 2423-2437.
#'
#' Davies, T.M., Jones, K. and Hazelton, M.L.
#' (2016), Symmetric adaptive smoothing regimens for estimation of the spatial
#' relative risk function, \emph{Computational Statistics & Data Analysis},
#' \bold{101}, 12-28.
#'
#' Diggle, P.J. (1985), A kernel method for smoothing
#' point process data, \emph{Journal of the Royal Statistical Society Series
#' C}, \bold{34}(2), 138-147.
#'
#' Hazelton, M.L. (2023), Shrinkage estimators of the spatial relative
#' risk function, \emph{Submitted for publication}.
#'
#' Hazelton, M.L. and Davies, T.M. (2009),
#' Inference based on kernel estimates of the relative risk function in
#' geographical epidemiology, \emph{Biometrical Journal}, \bold{51}(1),
#' 98-109.
#'
#' Kelsall, J.E. and Diggle, P.J. (1995a), Kernel estimation of
#' relative risk, \emph{Bernoulli}, \bold{1}, 3-16.
#'
#' Kelsall, J.E. and
#' Diggle, P.J. (1995b), Non-parametric estimation of spatial variation in
#' relative risk, \emph{Statistics in Medicine}, \bold{14}, 2335-2342.
#'
#' @examples
#' 
#' data(pbc)
#' pbccas <- split(pbc)$case
#' pbccon <- split(pbc)$control
#' h0 <- OS(pbc,nstar="geometric")
#' 
#' # Fixed (with tolerance contours)
#' pbcrr1 <- risk(pbccas,pbccon,h0=h0,tolerate=TRUE)
#' 
#' # Fixed shrinkage
#' pbcrr2 <- risk(pbccas,pbccon,h0=h0,shrink=TRUE,shrink.args=list(lambda=4))
#' 
#' # Asymmetric adaptive
#' pbcrr3 <- risk(pbccas,pbccon,h0=h0,adapt=TRUE,hp=c(OS(pbccas)/2,OS(pbccon)/2),
#'                tolerate=TRUE,davies.baddeley=0.05)
#' 
#' # Symmetric (pooled) adaptive
#' pbcrr4 <- risk(pbccas,pbccon,h0=h0,adapt=TRUE,tolerate=TRUE,hp=OS(pbc)/2,
#'                pilot.symmetry="pooled",davies.baddeley=0.05)
#' 
#' # Symmetric (case) adaptive; from two existing 'bivden' objects
#' f <- bivariate.density(pbccas,h0=h0,hp=2,adapt=TRUE,pilot.density=pbccas,
#'                        edge="diggle",davies.baddeley=0.05,verbose=FALSE) 
#' g <- bivariate.density(pbccon,h0=h0,hp=2,adapt=TRUE,pilot.density=pbccas,
#'                        edge="diggle",davies.baddeley=0.05,verbose=FALSE)
#' pbcrr5 <- risk(f,g,tolerate=TRUE,verbose=FALSE)
#' 
#' oldpar <- par(mfrow=c(2,2))
#' plot(pbcrr1,override.par=FALSE,main="Fixed")
#' plot(pbcrr2,override.par=FALSE,main="Fixed shrinkage")
#' plot(pbcrr3,override.par=FALSE,main="Asymmetric adaptive")
#' plot(pbcrr4,override.par=FALSE,main="Symmetric (pooled) adaptive") 
#' par(oldpar)
#' 
#' @export
risk <- function (f, g = NULL, log = TRUE, h0 = NULL, hp = h0, adapt = FALSE, shrink = FALSE, shrink.args = list(rescale = TRUE, type = c("lasso", "Bithell"), lambda = NA),
                   tolerate = FALSE, doplot = FALSE, pilot.symmetry = c("none", "f", "g", "pooled"), epsilon = 0, verbose = TRUE, ...) 
{
  if (is.null(g)) {
    if (!inherits(f, "ppp")) 
      stop("'f' must be an object of class 'ppp' if 'g' unsupplied")
    fm <- marks(f)
    if (!is.factor(fm)) 
      marks(f) <- fm <- factor(fm)
    if (nlevels(fm) != 2) 
      stop("'f' marks must be dichotomous if 'g' unsupplied")
    fs <- split(f)
    f <- fs[[1]]
    g <- fs[[2]]
  }
  else {
    fc <- class(f)
    gc <- class(g)
    if (!all(fc == gc)) 
      stop("'f' and 'g' must be of identical class")
    if (!(inherits(f, "ppp") || inherits(f, "bivden"))) 
      stop("'f' and 'g' must be of class 'ppp' or 'bivden'")
  }
  
  
  epsi <- epsilon[1]
  if (epsi > 0)
    warning("use of non-zero epsilon parameter is not recommended; use option 'shrink = TRUE' to employ either lasso or Bithell-type shrinkage")
  
  if (shrink) {
    epsi <- 0
    lambda <- shrink.args$lambda
    shrink.rescale <- shrink.args$rescale
    shrink.type <- shrink.args$type[1]
    if(is.null(lambda)) lambda <- NA
    if(is.null(shrink.rescale)) shrink.rescale <- TRUE
    if(is.null(shrink.type)) shrink.type <- "lasso"
    if(adapt) stop("shrinkage only implemented for fixed bandwidth estimators")
    if(is.numeric(lambda)&&(lambda<0)) stop("shrinkage parameter 'lambda' must be non-negative")
  }
  
  if (epsi < 0) 
    stop("invalid 'epsilon'; must be scalar and non-negative")
  
  
  if (inherits(f, "ppp")) {
    if (!identical_windows(Window(f), Window(g))) 
      stop("study windows for 'f' and 'g' must be identical")
    marks(f) <- NULL
    marks(g) <- NULL
    pooled <- suppressWarnings(superimpose(f, g))
    if (is.null(h0)) 
      h0 <- OS(pooled, nstar = sqrt(f$n * g$n))
    if (length(h0) == 1) {
      h0f <- h0g <- checkit(h0[1], "'h0[1]'")
    }
    else {
      h0f <- checkit(h0[1], "'h0[1]'")
      h0g <- checkit(h0[2], "'h0[2]'")
    }
    if (!adapt) {
      if (verbose) 
        message("Estimating case and control densities...", 
                appendLF = FALSE)
      
      fd <- bivariate.density(f, h0 = h0f, adapt = FALSE, ...)
      gd <- bivariate.density(g, h0 = h0g, adapt = FALSE, ...)
      
      if(shrink){
        if(h0f != h0g) stop("common case-control bandwidth required when employing shrinkage")
        
        h <- c(h0f,h0g)
        
        X1 <- fd$pp
        X2 <- gd$pp
        n1 <- npoints(X1)
        n2 <- npoints(X2)
        
        if (shrink.type=="lasso"){
          if (is.na(lambda)) lambda <- cv.RelRisk(X1,X2,h=h)$lambda
          case1 <- (fd$z-lambda/(n1*2*pi*h[1]^2))/(gd$z+lambda/(n2*2*pi*h[2]^2))
          case2 <- (fd$z+lambda/(n1*2*pi*h[1]^2))/(gd$z-lambda/(n2*2*pi*h[2]^2))
          logRR <- case1*0
          logRR[case1 > 1] <- log(case1[case1 > 1])
          logRR[1/case2 > 1] <- log(case2[1/case2 > 1])
          if (shrink.rescale) logRR <- logRR - log(spatstat.univar::integral(exp(logRR)*gd$z)) #+log(n2)
        } else if (shrink.type=="Bithell"){
          if (is.na(lambda)) lambda <- cv.RelRisk.Bithell(X1,X2,h=h)$lambda
          logRR <- log(fd$z+lambda/(n1*2*pi*h[1]^2)) - log(gd$z + lambda/(n1*2*pi*h[1]^2))
        } else {
          stop("invalid shrinkage type")
        }
      }
      
      if (verbose) 
        message("Done.")
      
    } # adaptive estimator here, still ppp
    else {
      if (is.null(hp)) 
        hp <- c(h0f, h0g)
      if (length(hp) == 1) {
        hfp <- hgp <- checkit(hp[1], "'hp[1]'")
      }
      else {
        hfp <- checkit(hp[1], "'hp[1]'")
        hgp <- checkit(hp[2], "'hp[2]'")
      }
      pilot.symmetry <- pilot.symmetry[1]
      pdat <- list()
      if (pilot.symmetry == "none") {
        pdat[[1]] <- f
        pdat[[2]] <- g
      }
      else if (pilot.symmetry == "f") {
        pdat[[1]] <- pdat[[2]] <- f
      }
      else if (pilot.symmetry == "g") {
        pdat[[1]] <- pdat[[2]] <- g
      }
      else if (pilot.symmetry == "pooled") {
        pdat[[1]] <- pdat[[2]] <- pooled
      }
      else {
        stop("invalid 'pilot.symmetry' argument")
      }
      if (verbose) 
        message("Estimating case density...", appendLF = FALSE)
      fd <- bivariate.density(f, h0 = h0f, hp = hfp, adapt = TRUE, 
                              pilot.density = pdat[[1]], verbose = FALSE, ...)
      if (verbose) 
        message("Done.\nEstimating control density...", 
                appendLF = FALSE)
      gd <- bivariate.density(g, h0 = h0g, hp = hgp, adapt = TRUE, 
                              pilot.density = pdat[[2]], verbose = FALSE, ...)
      if (verbose) 
        message("Done.")
    }
  }
  else { ##here bivden
    if (!compatible(f$z, g$z)) 
      stop("incompatible images in 'f' and 'g'... kernel estimates must be evaluated on identical domains")
    fd <- f
    gd <- g
    fda <- is.na(fd$gamma) || is.na(fd$geometric)
    gda <- is.na(gd$gamma) || is.na(gd$geometric)
    adapt <- switch(as.character(fda + gda), `0` = TRUE, 
                    `2` = FALSE, NA)
    if (is.na(adapt)) 
      stop("'f' and 'g' smoothed differently... must both be either fixed or adaptive")
    
    if(shrink){
      if(!is.na(f$gamma)||!is.na(g$gamma)) stop("shrinkage only implemented for fixed bandwidth estimates")
      h0f <- fd$h0
      h0g <- gd$h0
      if(f$h0!=g$h0) stop("common case-control bandwidth required when employing shrinkage")
      
      h <- c(h0f,h0g)
      X1 <- fd$pp
      X2 <- gd$pp
      n1 <- npoints(X1)
      n2 <- npoints(X2)
      
      if (shrink.type=="lasso"){
        if (is.na(lambda)) lambda <- cv.RelRisk(X1,X2,h=h)$lambda
        case1 <- (fd$z-lambda/(n1*2*pi*h[1]^2))/(gd$z+lambda/(n2*2*pi*h[2]^2))
        case2 <- (fd$z+lambda/(n1*2*pi*h[1]^2))/(gd$z-lambda/(n2*2*pi*h[2]^2))
        logRR <- case1*0
        logRR[case1 > 1] <- log(case1[case1 > 1])
        logRR[1/case2 > 1] <- log(case2[1/case2 > 1])
        if (shrink.rescale) logRR <- logRR - log(spatstat.univar::integral(exp(logRR)*gd$z))#+log(n2)
      } else if (shrink.type=="Bithell"){
        if (is.na(lambda)) lambda <- cv.RelRisk.Bithell(X1,X2,h=h)$lambda
        logRR <- log(fd$z+lambda/(n1*2*pi*h[1]^2)) - log(gd$z + lambda/(n1*2*pi*h[1]^2))
      } else {
        stop("invalid shrinkage type")
      }
    }
  }
  
  ## here final risk
  if (!shrink){
    eg <- epsi * max(gd$z)
    if (log) suppressWarnings(rr <- log(fd$z + eg) - log(gd$z + eg))
    else rr <- (fd$z + eg)/(gd$z + eg)
  } else {
    if (log) rr <- logRR
    else rr <- exp(logRR)
    #       fd <- list(z=f1,h0=h0,hp=NA,h=NA,him=NA,q=NA,gamma=NA,geometric=NA,pp=X1)
    # 	    gd <- list(z=f2,h0=h0,hp=NA,h=NA,him=NA,q=NA,gamma=NA,geometric=NA,pp=X2)
    # 	    class(fd) <- class(gd) <- "bivden"
  }
  
  
  if(tolerate & shrink) warning("tolerance contours not computed when 'shrink=TRUE'")
  
  ps <- NULL
  if (tolerate & !shrink) {
    if (verbose) 
      message("Calculating tolerance contours...", 
              appendLF = FALSE)
    if (adapt) 
      ps <- tol.asy.ada(fd, gd, 0.025, verbose = FALSE)$p
    else ps <- tol.asy.fix(fd, gd, gd, verbose = FALSE)$p
    if (verbose) 
      message("Done.")
  }
  if (doplot) {
    plot.im(rr, main = "", box = FALSE, ribargs = list(box = TRUE))
    axis(1)
    axis(2)
    box(bty = "l")
    plot(Window(fd$pp), add = TRUE)
    if (!is.null(ps)) 
      contour(fd$z$xcol, fd$z$yrow, t(as.matrix(ps)), levels = 0.05, 
              add = TRUE)
    return(invisible(NULL))
  }
  result <- list(rr = rr, f = fd, g = gd, P = ps)
  class(result) <- "rrs"
  return(result)
}
