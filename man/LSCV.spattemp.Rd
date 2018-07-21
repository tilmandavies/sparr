\name{LIK.spattemp}
\alias{LSCV.spattemp}
\alias{LIK.spattemp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation bandwidths for spatiotemporal kernel density estimates
}
\description{
Bandwidth selection for standalone spatiotemporal density/intensity based on either unbiased least squares cross-validation (LSCV) or likelihood (LIK) cross-validation, providing an isotropic scalar spatial bandwidth and a scalar temporal bandwidth.
}
\usage{
LIK.spattemp(pp, tt = NULL, tlim = NULL, sedge = c("uniform", "none"),
  tedge = sedge, parallelise = NA, start = NULL, verbose = TRUE)

LSCV.spattemp(pp, tt = NULL, tlim = NULL, sedge = c("uniform", "none"),
  tedge = sedge, sres = 64, tres = sres, parallelise = NA,
  start = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pp}{
An object of class \code{\link[spatstat]{ppp}} giving the spatial coordinates of the observations to be smoothed. Possibly marked with the time of each event; see argument \code{tt}.
}
  \item{tt}{
A numeric vector of equal length to the number of points in \code{pp}, giving the time corresponding to each spatial observation. If unsupplied, the function attempts to use the values in the \code{\link[spatstat]{marks}} attribute of the \code{\link[spatstat]{ppp.object}} in \code{pp}.
}
  \item{tlim}{
A numeric vector of length 2 giving the limits of the temporal domain over which to smooth. If supplied, all times in \code{tt} must fall within this interval (equality with limits allowed). If unsupplied, the function simply uses the range of the observed temporal values.
}
  \item{sedge}{
Character string dictating spatial edge correction. \code{"uniform"} (default) corrects based on evaluation grid coordinate. Setting \code{sedge="none"} requests no edge correction.
}
  \item{tedge}{
As \code{sedge}, for temporal edge correction.
}
  \item{sres}{
Numeric value > 0. Resolution of the [\code{sres} \eqn{\times}{x} \code{sres}] evaluation grid in the spatial margin.
}
  \item{tres}{
Numeric value > 0. Resolution of the evaluation points in the temporal margin as defined by the \code{tlim} interval. If unsupplied, the density is evaluated at integer values between \code{tlim[1]} and \code{tlim[2]}.
}
  \item{parallelise}{
Optional numeric argument to invoke parallel processing, by giving the number of CPU cores to use optimisation. This is only useful for larger data sets of many thousand observations. Experimental. Test your system first using \code{parallel::detectCores()} to identify the number of cores available to you.
}
  \item{start}{
Optional positive numeric vector of length 2 giving starting values for the internal call to \code{\link[stats]{optim}}, in the order of (<spatial bandwidth>, <temporal bandwidth>).
}
  \item{verbose}{
Logical value indicating whether to print a function progress bar to the console during evaluation.
}
}

\value{
A numeric vector of length 2 giving the jointly optimised spatial and temporal bandwidths (named \code{h} and \code{lambda} respectively).
}

\references{
Silverman, B.W. (1986), \emph{Density Estimation for Statistics and Data Analysis}, Chapman & Hall, New York.
}

\author{
T. M. Davies
}

\section{Warning}{
Leave-one-out CV for bandwidth selection in kernel
density estimation is notoriously unstable in practice and has a tendency to
produce rather small bandwidths in the fixed bandwidth case. Satisfactory bandwidths are not guaranteed
for every application. This method can also be computationally expensive for
large data sets and fine evaluation grid resolutions.
}

\seealso{
\code{\link{BOOT.spattemp}}, \code{\link{spattemp.density}}
}

\examples{
\donttest{

data(burk) # Burkitt's Uganda lymphoma data
burkcas <- burk$cases

hlam1 <- LSCV.spattemp(burkcas) #~9 secs
hlam2 <- LSCV.spattemp(burkcas,tlim=c(400,5800)) 
hlam3 <- LSCV.spattemp(burkcas,start=c(7,400))
rbind(hlam1,hlam2,hlam3)

hlam1 <- LIK.spattemp(burkcas) #~3 secs
hlam2 <- LIK.spattemp(burkcas,tlim=c(400,5800)) 
hlam3 <- LIK.spattemp(burkcas,start=c(7,400))
rbind(hlam1,hlam2,hlam3)

}
}
