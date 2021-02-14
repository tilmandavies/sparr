\name{spattemp.risk}
\alias{spattemp.risk}
\alias{rrst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatiotemporal relative risk/density ratio
}
\description{
Produces a spatiotemporal relative risk surface based on the ratio of two kernel estimates of spatiotemporal densities.
}
\usage{
spattemp.risk(f, g, log = TRUE, tolerate = FALSE, finiteness = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
An object of class \code{\link{stden}} representing the `case' (numerator) density estimate.
}
  \item{g}{
Either an object of class \code{\link{stden}}, or an object of class \code{\link{bivden}} for the `control' (denominator) density estimate. This object \bold{must} match the spatial (and temporal, if \code{\link{stden}}) domain of \code{f} completely; see `Details'.
}
  \item{log}{
Logical value indicating whether to return the log relative risk (default) or the raw ratio.
}
  \item{tolerate}{
Logical value indicating whether to compute and return asymptotic \eqn{p}-value surfaces for elevated risk; see `Details'.
}
  \item{finiteness}{
Logical value indicating whether to internally correct infinite risk (on the log-scale) to the nearest finite value to avoid numerical problems. A small extra computational cost is required.
}
  \item{verbose}{
Logical value indicating whether to print function progress during execution.
}
}
\details{
Fernando & Hazelton (2014) generalise the spatial relative risk function (e.g. Kelsall & Diggle, 1995) to the spatiotemporal domain. This is the implementation of their work, yielding the generalised log-relative risk function for \eqn{x\in W\subset R^2} and \eqn{t\in T\subset R}. It produces
\deqn{\hat{\rho}(x,t)=\log(\hat{f}(x,t))-\log(\hat{g}(x,t)),}
where \eqn{\hat{f}(x,t)} is a fixed-bandwidth kernel estimate of the spatiotemporal density of the cases (argument \code{f}) and \eqn{\hat{g}(x,t)} is the same for the controls (argument \code{g}).
\itemize{
\item When argument \code{g} is an object of class \code{\link{stden}} arising from a call to \code{\link{spattemp.density}}, the resolution, spatial domain, and temporal domain of this spatiotemporal estimate must match that of \code{f} exactly, else an error will be thrown.
\item When argument \code{g} is an object of class \code{\link{bivden}} arising from a call to \code{\link{bivariate.density}}, it is assumed the `at-risk' control density is static over time. In this instance, the above equation for the relative risk becomes \eqn{\hat{\rho}=\log(\hat{f}(x,t))+\log|T|-\log(g(x))}. The spatial density estimate in \code{g} must match the spatial domain of \code{f} exactly, else an error will be thrown.
\item The estimate \eqn{\hat{\rho}(x,t)} represents the joint or unconditional spatiotemporal relative risk over \eqn{W\times T}. This means that the raw relative risk \eqn{\hat{r}(x,t)=\exp{\hat{\rho}(x,t)}} integrates to 1 with respect to the control density over space and time: \eqn{\int_W \int_T r(x,t)g(x,t) dt dx = 1}. This function also computes the \bold{conditional} spatiotemporal relative risk at each time point, namely
\deqn{\hat{\rho}(x|t)=\log{\hat{f}(x|t)}-\log{\hat{g}(x|t)},}
where \eqn{\hat{f}(x|t)} and \eqn{\hat{g}(x|t)} are the conditional densities over space of the cases and controls given a specific time point \eqn{t} (see the documentation for \code{\link{spattemp.density}}). In terms of normalisation, we therefore have \eqn{\int_W r(x|t)g(x|t) dx = 1}. In the case where \eqn{\hat{g}} is static over time, one may simply replace \eqn{\hat{g}(x|t)} with \eqn{\hat{g}(x)} in the above.
\item Based on the asymptotic properties of the estimator, Fernando & Hazelton (2014) also define the calculation of tolerance contours for detecting statistically significant fluctuations in such spatiotemporal log-relative risk surfaces. This function can produce the required \eqn{p}-value surfaces by setting \code{tolerate = TRUE}; and if so, results are returned for both the unconditional (x,t) and conditional (x|t) surfaces. See the examples in the documentation for \code{\link{plot.rrst}} for details on how one may superimpose contours at specific \eqn{p}-values for given evaluation times \eqn{t} on a plot of relative risk on the spatial margin. %Calculation of these \eqn{p}-values requires a common estimate of the space-time heterogeneity, which is performed internally by pooling the case and control data, and making use of the geometric means of the spatial- and temporal-margin bandwidths in \code{f} and \code{g} (typically, these are the same for both case and control densities anyway).
}
}

\value{
An object of class ``\code{rrst}''. This is effectively a list with the following members:
\item{rr}{
A named (by time-point) list of pixel \code{\link[spatstat.geom]{im}}ages corresponding to the joint spatiotemporal relative risk over space at each discretised time.
}
\item{rr.cond}{
A named list of pixel \code{\link[spatstat.geom]{im}}ages corresponding to the conditional spatial relative risk given each discretised time.
}
\item{P}{
A named list of pixel \code{\link[spatstat.geom]{im}}ages of the \eqn{p}-value surfaces testing for elevated risk for the joint estimate. If \code{tolerate = FALSE}, this will be \code{NULL}.
}
\item{P.cond}{
As above, for the conditional relative risk surfaces.
}
\item{f}{
A copy of the object \code{f} used in the initial call.
}
\item{g}{
As above, for \code{g}.
}
\item{tlim}{
A numeric vector of length two giving the temporal bound of the density estimate.
}
}

\references{
Fernando, W.T.P.S. and Hazelton, M.L. (2014), Generalizing the spatial relative risk function, \emph{Spatial and Spatio-temporal Epidemiology}, \bold{8}, 1-10.\cr\cr
}
\author{
T.M. Davies
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{spattemp.density}}, \code{\link{spattemp.slice}}, \code{\link{bivariate.density}}
}
\examples{
\donttest{
data(fmd)
fmdcas <- fmd$cases
fmdcon <- fmd$controls

f <- spattemp.density(fmdcas,h=6,lambda=8) # stden object as time-varying case density
g <- bivariate.density(fmdcon,h0=6) # bivden object as time-static control density
rho <- spattemp.risk(f,g,tolerate=TRUE) 
print(rho)

par(mfrow=c(2,3))
plot(rho$f$spatial.z,main="Spatial margin (cases)") # spatial margin of cases
plot(rho$f$temporal.z,main="Temporal margin (cases)") # temporal margin of cases
plot(rho$g$z,main="Spatial margin (controls)") # spatial margin of controls
plot(rho,tselect=50,type="conditional",tol.args=list(levels=c(0.05,0.0001),
     lty=2:1,lwd=1:2),override.par=FALSE)
plot(rho,tselect=100,type="conditional",tol.args=list(levels=c(0.05,0.0001),
     lty=2:1,lwd=1:2),override.par=FALSE)
plot(rho,tselect=200,type="conditional",tol.args=list(levels=c(0.05,0.0001),
     lty=2:1,lwd=1:2),override.par=FALSE)
}
}
