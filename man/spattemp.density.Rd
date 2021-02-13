\name{spattemp.density}
\alias{spattemp.density}
\alias{stden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatiotemporal kernel density estimation
}
\description{
Provides a fixed-bandwidth kernel estimate of continuous spatiotemporal data.
}
\usage{
spattemp.density(pp, h = NULL, tt = NULL, lambda = NULL,
  tlim = NULL, sedge = c("uniform", "none"), tedge = sedge, 
  sres = 128, tres = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pp}{
An object of class \code{\link[spatstat.geom]{ppp}} giving the spatial coordinates of the observations to be smoothed. Possibly marked with the time of each event; see argument \code{tt}.
}
  \item{h}{
Fixed bandwidth to smooth the spatial margin. A numeric value > 0. If unsupplied, the oversmoothing bandwidth is used as per \code{\link{OS}}.
}
  \item{tt}{
A numeric vector of equal length to the number of points in \code{pp}, giving the time corresponding to each spatial observation. If unsupplied, the function attempts to use the values in the \code{\link[spatstat.geom]{marks}} attribute of the \code{\link[spatstat.geom:ppp]{ppp.object}} in \code{pp}.
}
  \item{lambda}{
Fixed bandwidth to smooth the temporal margin; a numeric value > 0. If unsupplied, the function internally computes the Sheather-Jones bandwith using \code{\link[stats]{bw.SJ}} (Sheather & Jones, 1991).
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
  \item{verbose}{
Logical value indicating whether to print a function progress bar to the console during evaluation.
}
}

\details{
This function produces a fixed-bandwidth kernel estimate of a single spatiotemporal density, with isotropic smoothing in the spatial margin, as per Fernando & Hazelton (2014). Estimates may be edge-corrected for an irregular spatial study window \emph{and} for the bounds on the temporal margin as per \code{tlim}; this edge-correction is performed in precisely the same way as the \code{"uniform"} option in \code{\link{bivariate.density}}.

Specifically, for \eqn{n} trivariate points in space-time (\code{pp}, \code{tt}, \code{tlim}), we have
\deqn{\hat{f}(x,t)=n^{-1}\sum_{i=1}^{n}h^{-2}\lambda^{-1}K((x-x_i)/h)L((t-t_i)/\lambda)/(q(x)q(t)),}
where \eqn{x\in W\subset R^2} and \eqn{t\in T\subset R}; \eqn{K} and \eqn{L} are the 2D and 1D Gaussian kernels controlled by fixed bandwidths \eqn{h} (\code{h}) and \eqn{\lambda} (\code{lambda}) respectively; and \eqn{q(x)=\int_W h^{-2}K((u-x)/h)du} and \eqn{q(t)=\int_T \lambda^{-1}L((w-t)/\lambda)dw} are optional edge-correction factors (\code{sedge} and \code{tedge}).

The above equation provides the \emph{joint} or \emph{unconditional} density at a given space-time location \eqn{(x,t)}. In addition to this, the function also yields the \emph{conditional} density at each grid time, defined as
\deqn{\hat{f}(x|t)=\hat{f}(x,t)/\hat{f}(t),}
where \eqn{\hat{f}(t)=n^{-1}\sum_{i=1}^{n}\lambda^{-1}L((t-t_i)/\lambda)/q(t)} is the univariate kernel estimate of the temporal margin. Normalisation of the two versions \eqn{\hat{f}(x,t)} and \eqn{\hat{f}(x|t)} is the only way they differ. Where in the unconditional setting we have \eqn{\int_W\int_T\hat{f}(x,t)dt dx=1}, in the conditional setting we have \eqn{\int_W\hat{f}(x|t) dx=1} for all \eqn{t}. See Fernando & Hazelton (2014) for further details and practical reasons as to why we might prefer one over the other in certain situations.

The objects returned by this function (see `Value' below) are necessary for kernel estimation of spatiotemporal relative risk surfaces, which is performed by \code{\link{spattemp.risk}}.
}


\value{
An object of class \code{"stden"}. This is effectively a list with the following components:
\item{z}{
A named (by time-point) list of pixel \code{\link[spatstat.geom]{im}}ages corresponding to the joint spatiotemporal density over space at each discretised time.
}
\item{z.cond}{
A named (by time-point) list of pixel \code{\link[spatstat.geom]{im}}ages corresponding to the conditional spatial density given each discretised time.
}
\item{h}{
The scalar bandwidth used for spatial smoothing.
}
\item{lambda}{
The scalar bandwidth used for temporal smoothing.
}
\item{tlim}{
A numeric vector of length two giving the temporal bound of the density estimate.
}
\item{spatial.z}{
A pixel \code{\link[spatstat.geom]{im}}age giving the overall spatial margin as a single 2D density estimate (i.e. ignoring time).
}
\item{temporal.z}{
An object of class \code{\link[stats]{density}} giving the overall temporal margin as a single 1D density estimate (i.e. ignoring space).
}
\item{qs}{
A pixel \code{\link[spatstat.geom]{im}}age giving the edge-correction surface for the spatial margin. \code{NULL} if \code{sedge = "none"}.
}
\item{qt}{
A numeric vector giving the edge-correction weights for the temporal margin. \code{NULL} if \code{tedge = "none"}.
}
\item{pp}{
A \code{\link[spatstat.geom:ppp]{ppp.object}} of the spatial data passed to the argument of the same name in the initial function call, with \code{\link[spatstat.geom]{marks}} of the observation times.
}
\item{tgrid}{
A numeric vector giving the discretised time grid at which the spatiotemporal density was evaluated (matches the names of \code{z} and \code{z.cond}).
}
}

\references{
Duong, T. (2007), ks: Kernel Density Estimation and Kernel Discriminant Analysis for Multivariate Data in R, \emph{Journal of Statistical Software}, \bold{21}(7), 1-16.\cr\cr
Fernando, W.T.P.S. and Hazelton, M.L. (2014), Generalizing the spatial relative risk function, \emph{Spatial and Spatio-temporal Epidemiology}, \bold{8}, 1-10.\cr\cr
Kelsall, J.E. and Diggle, P.J. (1995), Kernel estimation of relative risk, \emph{Bernoulli}, \bold{1}, 3-16.\cr\cr
Sheather, S. J. and Jones, M. C. (1991), A reliable data-based bandwidth selection method for kernel density estimation. Journal of the Royal Statistical Society Series B, \bold{53}, 683-690.\cr\cr
Silverman, B.W. (1986), \emph{Density Estimation for Statistics and Data Analysis}, Chapman & Hall, New York.
}

\author{
T.M. Davies
}

\seealso{
\code{\link{bivariate.density}}, \code{\link{spattemp.risk}}, \code{\link{spattemp.slice}}
}

\examples{
data(burk)
burkcas <- burk$cases

burkden1 <- spattemp.density(burkcas,tres=128)
summary(burkden1)

\donttest{
hlam <- LIK.spattemp(burkcas,tlim=c(400,5900),verbose=FALSE)
burkden2 <- spattemp.density(burkcas,h=hlam[1],lambda=hlam[2],tlim=c(400,5900),tres=256)
tims <- c(1000,2000,3500)
par(mfcol=c(2,3))
for(i in tims){ 
  plot(burkden2,i,override.par=FALSE,fix.range=TRUE,main=paste("joint",i))
  plot(burkden2,i,"conditional",override.par=FALSE,main=paste("cond.",i))
}
}
}

%Evaluation of the spatiotemporal density itself is direct, limiting effective support of the Gaussian kernel to plus or minus 4 bandwidths, based on the implementation of 3D density estimation in the \code{\link[ks]{ks}} package (Duong, 2007). 
