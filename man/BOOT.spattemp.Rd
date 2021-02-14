\name{BOOT.spattemp}
\alias{BOOT.spattemp}

\title{
Bootstrap bandwidths for a spatiotemporal kernel density estimate
}
\description{
Bandwidth selection for standalone spatiotemporal density/intensity based on bootstrap estimation of the MISE, providing an isotropic scalar spatial bandwidth and a scalar temporal bandwidth.
}
\usage{
BOOT.spattemp(pp, tt = NULL, tlim = NULL, eta = NULL, nu = NULL,
  sedge = c("uniform", "none"), tedge = sedge, ref.density = NULL,
  sres = 64, tres = sres, start = NULL, verbose = TRUE)
}

\arguments{
  \item{pp}{
An object of class \code{\link[spatstat.geom]{ppp}} giving the spatial coordinates of the observations to be smoothed. Possibly marked with the time of each event; see argument \code{tt}.
}
  \item{tt}{
A numeric vector of equal length to the number of points in \code{pp}, giving the time corresponding to each spatial observation. If unsupplied, the function attempts to use the values in the \code{\link[spatstat.geom]{marks}} attribute of the \code{\link[spatstat.geom:ppp]{ppp.object}} in \code{pp}.
}
  \item{tlim}{
A numeric vector of length 2 giving the limits of the temporal domain over which to smooth. If supplied, all times in \code{tt} must fall within this interval (equality with limits allowed). If unsupplied, the function simply uses the range of the observed temporal values.
}
  \item{eta}{
Fixed scalar bandwidth to use for the spatial margin of the reference density estimate; if \code{NULL} it is calculated as the oversmoothing bandwidth of \code{pp} using \code{\link{OS}}. Ignored if \code{ref.density} is supplied. See `Details'.
}
  \item{nu}{
Fixed scalar bandwidth to use for the temporal margin of the reference density estimate; if \code{NULL} it is calculated from \code{tt} using the univariate version of Terrell's (1990) oversmoothing principle. Ignored if \code{ref.density} is supplied. See `Details'.
}
  \item{sedge}{
Character string dictating spatial edge correction. \code{"uniform"} (default) corrects based on evaluation grid coordinate. Setting \code{sedge="none"} requests no edge correction.
}
  \item{tedge}{
As \code{sedge}, for temporal edge correction.
}
  \item{ref.density}{
Optional. An object of class \code{\link{stden}} giving the reference density from which data is assumed to originate in the bootstrap. Must be spatially edge-corrected if \code{sedge = "uniform"}.
}
  \item{sres}{
Numeric value > 0. Resolution of the [\code{sres} \eqn{\times}{x} \code{sres}] evaluation grid in the spatial margin.
}
  \item{tres}{
Numeric value > 0. Resolution of the evaluation points in the temporal margin as defined by the \code{tlim} interval. If unsupplied, the density is evaluated at integer values between \code{tlim[1]} and \code{tlim[2]}.
}
  \item{start}{
Optional positive numeric vector of length 2 giving starting values for the internal call to \code{\link[stats]{optim}}, in the order of (<spatial bandwidth>, <temporal bandwidth>).
}
  \item{verbose}{
Logical value indicating whether to print a function progress bar to the console during evaluation.
}
}

\details{
For a spatiotemporal kernel density estimate \eqn{\hat{f}} defined on \eqn{W x T \in R^3}, the mean integrated squared error (MISE) is given by \eqn{E[\int_W \int_T (\hat{f}(x,t) - f(x,t))^2 dt dx]}, where \eqn{f} is the corresponding true density. Given observed spatiotemporal locations \eqn{X} (arguments \code{pp} and \code{tt}) of \eqn{n} observations, this function finds the scalar spatial bandwidth \eqn{h} and scalar temporal bandwidth \eqn{\lambda} that jointly minimise
\deqn{E^*[\int_W \int_T (\hat{f}^*(x,t) - \hat{f}(x,t))^2 dt dx],}
where \eqn{\hat{f}(x,t)} is a density estimate of \eqn{X} constructed with `reference' bandwidths \eqn{\eta} (spatial; argument \code{eta}) and \eqn{\nu} (temporal; argument \code{nu}); \eqn{\hat{f}^*(x,t)} is a density estimate using bandwidths \eqn{h} and \eqn{\lambda} of \eqn{n} observations \eqn{X^*} generated from \eqn{\hat{f}(x,t)}. The notation \eqn{E^*} denotes expectation with respect to the distribution of the \eqn{X^*}. The user may optionally supply \code{ref.density} as an object of class \code{\link{stden}}, which must be evaluated on the same spatial and temporal domains \eqn{W} and \eqn{T} as the data (arguments \code{pp}, \code{tt}, and \code{tlim}). In this case, the reference bandwidths are extracted from this object, and \code{eta} and \code{nu} are ignored.

This function is based on an extension of the theory of Taylor (1989) to the spatiotemporal domain and to cope with the inclusion of edge-correction factors. No resampling is necessary due to the theoretical properties of the Gaussian kernel.
}

\value{
A numeric vector of length 2 giving the jointly optimised spatial and temporal bandwidths (named \code{h} and \code{lambda} respectively).
}

\references{
Taylor, C.C. (1989) Bootstrap choice of the smoothing parameter in kernel density estimation, \emph{Biometrika}, \bold{76}, 705-712.
}

\author{
T. M. Davies
}


\section{Warning}{
Bootstrapping for spatiotemporal bandwidth selection for spatiotemporal data is very computationally demanding. Keeping \code{verbose = TRUE} offers an indication of the computational burden by printing each pair of bandwidths at each iteration of the \code{\link{optim}}isation routine. The `Examples' section also offers some rough indications of evaluation times on this author's local machine.
}

\seealso{
\code{\link{LSCV.spattemp}}, \code{\link{spattemp.density}}
}

\examples{
\donttest{

data(burk) # Burkitt's Uganda lymphoma data
burkcas <- burk$cases

#~85 secs
hlam1 <- BOOT.spattemp(burkcas) 

#~75 secs. Widen time limits, reduce ref. bw.
hlam2 <- BOOT.spattemp(burkcas,tlim=c(400,5800),eta=8,nu=450) 

#~150 secs. Increase ref. bw., custom starting vals
hlam3 <- BOOT.spattemp(burkcas,eta=20,nu=800,start=c(7,400)) 

rbind(hlam1,hlam2,hlam3)
}
}
