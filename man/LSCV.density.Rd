\name{LSCV.density}
\alias{LSCV.density}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Leave-one-out least-squares cross-validation (LSCV) bandwidth selector
}
\description{
Isotropic fixed bandwidth selection for standalone 2D density/intensity based on classical unbiased cross-validation
}
\usage{
LSCV.density(pp, hlim = NULL, hseq = NULL, resolution = 64,
             edge = TRUE, auto.optim = TRUE, seqres = 30,
             parallelise = NULL, verbose = TRUE, type = "spatial",
             lambdalim = NULL, lambdaseq = NULL, tlim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pp}{
An object of class \code{\link[spatstat]{ppp}} giving the observed 2D data to be smoothed.
}
  \item{hlim}{
An optional vector of length 2 giving the limits of the optimisation routine with respect to the bandwidth. If unspecified, the function attempts to choose this automatically.
}
  \item{hseq}{
An optional increasing sequence of bandwidth values at which to manually evaluate the optimisation criterion. Used only in the case \code{(!auto.optim && is.null(hlim))}.  
}
  \item{resolution}{
Spatial grid size; the optimisation will be based on a [\code{resolution} \eqn{\times}{x} \code{resolution}] density estimate.
}
  \item{edge}{
Logical value indicating whether to edge-correct the density estimates used.  
}
  \item{auto.optim}{
Logical value indicating whether to automate the numerical optimisation using \code{\link{optimise}}. If \code{FALSE}, the optimisation criterion is evaluated over \code{hseq} (if supplied), or over a seqence of values controlled by \code{hlim} and \code{seqres}.  
}
  \item{seqres}{
Optional resolution of an increasing sequence of bandwidth values. Only used if \code{(!auto.optim && is.null(hseq))}.
}
  \item{parallelise}{
Numeric argument to invoke parallel processing, giving the number of CPU cores to use when \code{!auto.optim}. Experimental. Test your system first using \code{parallel::detectCores()} to identify the number of cores available to you.
}
  \item{verbose}{
Logical value indicating whether to provide function progress commentary.  
}
  \item{type}{
Unimplemented. Future argument for spatiotemporal bandwidth selection.  
}
  \item{lambdalim}{
Unimplemented. Future argument for spatiotemporal bandwidth selection.  
}
  \item{lambdaseq}{
Unimplemented. Future argument for spatiotemporal bandwidth selection.
}
  \item{tlim}{
Unimplemented. Future argument for spatiotemporal bandwidth selection.
}

}

\value{
A single numeric value of the estimated bandwidth (if \code{auto.optim = TRUE}). Otherwise, a list of two numeric vectors of equal length giving the bandwidth sequence (as \code{hs}) and corresponding CV function value (as \code{CV}).

}
\references{
Silverman, B.W. (1986), \emph{Density Estimation for Statistics and Data Analysis}, Chapman & Hall, New York. \cr\cr
Wand, M.P. and Jones, C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.
}
\author{
T. M. Davies
}

\section{Warning}{
Leave-one-out LSCV for bandwidth selection in kernel density estimation is notoriously unstable in practice and has a tendency to produce rather small bandwidths. Satisfactory bandwidths are not guaranteed for every application. This method can also be computationally expensive for large data sets and fine evaluation grid resolutions. The user may need to experiment with adjusting \code{hlim} to find a suitable minimum.
}

\seealso{
Functions for bandwidth selection in package \code{\link{spatstat}}: \code{\link[spatstat]{bw.diggle}}; \code{\link[spatstat]{bw.ppl}}; \code{\link[spatstat]{bw.scott}}; \code{\link[spatstat]{bw.frac}}.
}
\examples{
## To be filled
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
