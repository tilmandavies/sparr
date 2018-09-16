\name{BOOT.density}
\alias{BOOT.density}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap bandwidth for a spatial kernel density estimate
}
\description{
Isotropic fixed or global (for adaptive) bandwidth selection for a standalone 2D density based on bootstrap estimation of the MISE.
}
\usage{
BOOT.density(pp, hlim = NULL, eta = NULL, type = c("fixed", "adaptive"),
  hp = NULL, edge = c("uniform", "none"), ref.density = NULL,
  resolution = 64, rmdiag = TRUE, sim.adapt = list(N = 50, B = 100,
  dimz = 64, objective = FALSE), parallelise = NA, verbose = TRUE, ...)
}

\arguments{
  \item{pp}{
An object of class \code{\link[spatstat]{ppp}} giving the observed 2D data to be smoothed.
}
  \item{hlim}{
An optional vector of length 2 giving the limits of the optimisation routine with respect to the bandwidth. If \code{NULL}, the function attempts to choose this automatically.
}
  \item{eta}{
Fixed scalar bandwidth to use for the reference density estimate; if \code{NULL} it is calculated as the oversmoothing bandwidth of \code{pp} using \code{\link{OS}}. Ignored if \code{ref.density} is supplied. See `Details'. 
}
  \item{type}{
A character string indicating selection type. Either \code{"fixed"} (default) for selection of a constant bandwidth for the fixed-bandwidth estimator based on theory extended from results in Taylor (1989); or \code{"adaptive"} for selection of the global bandwidth for an adaptive kernel density. See `Details'.
}
  \item{hp}{
Pilot bandwidth used for adaptive estimates in the bootstrap; see the argument of the same tag in \code{\link{bivariate.density}}. Ignored when \code{type = "fixed"} or when \code{ref.density} is supplied.
}
  \item{edge}{
Character string dictating edge correction for the bootstrapped estimates. \code{"uniform"} (default) corrects based on evaluation grid coordinate. Setting \code{edge="none"} requests no edge correction.
}
  \item{ref.density}{
Optional. An object of class \code{\link{bivden}} giving the reference density from which data will be generated. Based on theory, this must be a fixed-bandwidth estimate if \code{type = "fixed"}; see `Details'. Must be edge-corrected if \code{edge = "uniform"}.
}
  \item{resolution}{
Spatial grid size; the optimisation will be based on a [\code{resolution} \eqn{\times}{x} \code{resolution}] density estimate.
}
  \item{rmdiag}{
Logical control value for removal of mirrored evaluation points as suggested by Taylor (1989) in the theoretical expression of the fixed-bandwidth MISE estimate. See `Details'. Ignored when \code{type = "adaptive"}
}
  \item{sim.adapt}{
List of control values for bootstrap simulation in the adaptive case; see `Details'. Ignored when \code{type = "fixed"}.
}
  \item{parallelise}{
Optional numeric argument to reduce computation time by invoking parallel processing, by giving the number of CPU cores to use in either evaluation (fixed) or in the actual bootstrap replicate generation (adaptive). Experimental. Test your system first using \code{parallel::detectCores()} to identify the number of cores available to you.
}
  \item{verbose}{
Logical value indicating whether to print function progress during execution.
}
  \item{\dots}{
Optional arguments controlling scaling to be passed to \code{\link{multiscale.density}} for the adaptive bootstrap; ignored when \code{type = "fixed"}.
}
}
\details{
For a 2D kernel density estimate \eqn{\hat{f}} defined on \eqn{W \in R^2}, the mean integrated squared error (MISE) is given by \eqn{E[\int_W (\hat{f}(x) - f(x))^2 dx]}, where \eqn{f} is the corresponding true density. Given an observed data set \eqn{X} (argument \code{pp}) of \eqn{n} observations, this function finds the bandwidth \eqn{h} that minimises
\deqn{E^*[\int_W (\hat{f}^*(x) - \hat{f}(x))^2 dx],}
where \eqn{\hat{f}(x)} is a density estimate of \eqn{X} constructed with `reference' bandwidth \eqn{\eta} (argument \code{eta} or \code{ref.density}), and \eqn{\hat{f}^*(x)} is a density estimate using bandwidth \eqn{h} of \eqn{n} observations \eqn{X^*} generated from \eqn{\hat{f}(x)}. The notation \eqn{E^*} denotes expectation with respect to the distribution of the \eqn{X^*}.

\describe{
\item{\bold{Fixed}}{
When \code{type = "fixed"}, the function assumes you want to select a constant bandwidth for use with the fixed-bandwith density estimator. This implementation is based on extending the remarkable results of Taylor (1989) (see also Sain et al., 1994), who demonstrates that when the Gaussian kernel is being used, we can find the optimal \eqn{h} with respect to the aforementioned bootstrap-estimated MISE without any actual resampling. This implementation extends these results to the bivariate setting, and allows for edge-correction of both the reference and bootstrap densities.
\itemize{
\item Taylor (1989) does not distinguish between the reference bandwidth \eqn{\eta} and the target of optimisation, \eqn{h}, thus allowing the reference bandwidth to vary alongside the target in the optimisation. This is not optimal, and this function always assumes a static reference bandwidth. Hall et al. (1992) indicate that a generous amount of smoothing is to be preferred in the reference density (hence the default \code{eta} set using \code{\link{OS}}).
\item If \code{ref.density} is supplied, it \bold{must} be a fixed-bandwidth density estimate as an object of class \code{\link{bivden}} for validity of the theory. Edge-correction must be present if \code{edge = "uniform"}; and it must be evaluated on the same spatial domain as dictated by \code{Window(pp)} and \code{resolution}. If unsupplied, the function internally computes an appropriate fixed-bandwidth density estimate using \code{eta} as the reference bandwidth.
\item Finally, Taylor (1989) argues it is preferable to avoid summation at identical evaluation grid points in the expression for the optimal bandwidth, which is performed when \code{rmdiag = TRUE}. Setting \code{rmdiag = FALSE} disables this correction.
}
}
\item{\bold{Adaptive}}{
When \code{type = "adaptive"}, the function assumes you want to select a global bandwidth (argument \code{h0} in \code{\link{bivariate.density}}) for use in 2D adaptive kernel density estimation.
\itemize{
\item An expression similar to Taylor (1989) is not possible for the adaptive estimator. Thus, in the adaptive setting, the optimal bootstrap bandwidth is calculated by brute force as was performed in Davies and Baddeley (2018) by taking advantage of the multiscale estimation theory implemented in \code{\link{multiscale.density}}. The value that minimises an interpolating cubic spline of the estimated MISE on bandwidth is identified as the optimal global bandwidth.
\item The user can pass either a fixed or adaptive \code{bivden} object to \code{ref.density}. If this is the case, \code{hp} is ignored and the pilot bandwidth for each iteration of the bootstrap in estimation of the \eqn{\hat{f}^*(x)} uses \code{ref.density$hp} (if \code{ref.density} is adaptive) or \code{ref.density$h0} (if \code{ref.density} is fixed). When \code{ref.density} is unsupplied, the function uses a fixed-bandwidth kernel estimate with bandwidth \code{eta} as the reference density, and if additionally \code{hp} is unsupplied, the same value \code{eta} is used for the constant pilot bandwidth.
\item Control over the bootstrap is achieved with four optional named arguments passed as a list to \code{sim.adapt}. \code{N} controls the number of bootstrap iterates per bandwidth; \code{B} controls the resolution of the sequence of bandwidths trialled (i.e. between \code{hlim[1]} and \code{hlim[2]}); \code{dimz} specifies the resolution of the bandwidth axis in the trivariate convolution evaluated by \code{\link{multiscale.density}}; and \code{objective} specifies whether to return the set of estimated MISEs for all bandwidths (nice to plot), or merely the optimal bandwidth (see `Value').
\item The \code{\dots} are intended for any relevant optional arguments to be passed to the internal call to \code{\link{multiscale.density}}, such as \code{gamma.scale} or \code{trim}.
}
}
}
}

\value{
The optimal fixed or global (for adaptive) scalar bandwidth. If \code{simargs$objective = TRUE} for the adaptive bootstrap, the return object is instead a \eqn{[}\code{simargs$B} \eqn{x 2]} matrix, with the first column giving the trialled bandwidth and the second giving the corresponding value of the estimated bootstrap MISE.
}

\references{
Davies, T.M. and Baddeley A. (2018), Fast computation of spatially adaptive kernel estimates, \emph{Statistics and Computing}, \bold{28}(4), 937-956.\cr\cr
Hall, P., Marron, J.S. and Park, B.U. (1992) Smoothed cross-validation, \emph{Probability Theory and Related Fields}, \bold{92}, 1-20.\cr\cr
Sain, S.R., Baggerly, K.A. and Scott, D.W. (1994) Cross-validation of multivariate densities, \emph{Journal of the American Statistical Association}, \bold{89}, 807-817.\cr\cr
Taylor, C.C. (1989) Bootstrap choice of the smoothing parameter in kernel density estimation, \emph{Biometrika}, \bold{76}, 705-712.
}

\author{
T.M. Davies
}

\section{Warning}{
Even with the implemented computational tricks, bootstrapping for bandwidth selection for spatial data is still computationally demanding, especially for adaptive kernel estimates. The user can reduce this time by keeping the evaluation grid at modest \code{resolution}s, and experimenting with parallelising the internal loops via \code{parallelise}. The `Examples' section offers some rough indications of evaluation times on this author's local machine.
}

\seealso{
\code{\link{bivariate.density}}, \code{\link{OS}}, \code{\link{multiscale.density}}
}

\examples{
\donttest{
data(pbc)

## Fixed bandwidth selection ##
BOOT.density(pbc) # ~20 secs
BOOT.density(pbc,eta=OS(pbc)/2) # halve default reference bandwidth
BOOT.density(pbc,eta=OS(pbc)*2) # double default reference bandwidth

# supplying pre-defined reference density as fixed-bandwidth 'bivden' object
pbcfix <- bivariate.density(pbc,h0=2.5,resolution=64)
system.time(hfix <- BOOT.density(pbc,ref.density=pbcfix,parallelise=4)) # parallelisation; 14 secs
hfix

## Global (for adaptive) bandwidth selection ##
# ~200 secs next line; use 'parallelise' for speedup
system.time(hada <- BOOT.density(pbc,type="adaptive")) # minimal usage for adaptive bootstrap
hada

# ~80 secs next line. Set custom h limits; increase reference bandwidth;
#    set custom pilot bandwidth; return objective function
system.time(hada <- BOOT.density(pbc,hlim=c(0.9,8),eta=3.5,type="adaptive",
                                 hp=OS(pbc)/2,parallelise=6,
                                 sim.adapt=list(objective=TRUE)))
hada[which.min(hada[,2]),1]
plot(hada);abline(v=hada[which.min(hada[,2]),1],col=2)
}
}
