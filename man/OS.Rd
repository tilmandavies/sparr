\name{OS}
\Rdversion{1.1}
\alias{OS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Oversmoothing (OS) bandwidth selector
}
\description{
Provides the isotropic fixed bandwidth for a 2D point pattern based on the maximal smoothing (oversmoothing) principle of Terrell (1990).
}
\usage{
OS(pp, nstar = c("npoints", "geometric"),
   scaler = c("silverman", "IQR", "sd", "var"))
}
\arguments{
  \item{pp}{
An object of class \code{\link[spatstat]{ppp}} giving the observed 2D data to be smoothed.
}
  \item{nstar}{
Optional. Controls the value to use in place of the number of observations \emph{n} in the oversmoothing formula. Either a character string, \code{"npoints"} (default) or \code{"geometric"}, or a positive numeric value. See `Details'. 
}
  \item{scaler}{
Optional. Controls the value for a scalar representation of the 2D scale of the data. Either a character string, \code{"silverman"} (default), \code{"IQR"}, \code{"sd"}, or \code{"var"}; or a positive numeric value. See `Details'.}
}

\details{
This function calculates a scalar smoothing bandwidth for kernel density estimates of 2-dimensional data: the ``maximal amount of smoothing compatible with the estimated scale of the observed data''. See Terrell (1990).
\describe{
\item{\bold{Effective sample size}}{
The formula requires a sample size, and this can be minimally tailored via \code{nstar}. By default, the function simply uses the number of observations in \code{pp}: \code{nstar = "npoints"}. If \code{pp} is a \code{\link[spatstat]{ppp.object}} with factor-valued \code{\link[spatstat]{marks}}, then the user has the option of using \code{nstar = "geometric"}, which sets the sample size used in the formula to the geometric mean of the counts of observations of each mark. This can be useful for e.g. relative risk calculations, see Davies and Hazelton (2010). Alternatively, the user can specify their own value by simply supplying a single positive numeric value to \code{nstar}.
}
\item{\bold{Spatial scale}}{
For isotropic smoothing, we require a scalar estimate of the spatial scale of the data. One may use the `robust' estimate of standard deviation found by a weighted mean of the interquartile ranges of the \eqn{x}- and \eqn{y}-coordinates of the data respectively (\code{scaler = "IQR"}). Two other options are the raw mean of the coordinate-wise standard deviations (\code{scaler = "sd"}), or the square root of the mean of the two variances (\code{scaler = "var"}). A fourth option, \code{scaler = "silverman"} (default), sets the scaling constant to be the minimum of the \code{"IQR"} and \code{"sd"} options; see Silverman (1986), p. 47. Alternatively, like \code{nstar}, the user can specify their own value by simply supplying a single positive numeric value to \code{scaler}.
}
}
}

\value{
A single numeric value of the estimated bandwidth.
}

\references{
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.\cr\cr
Terrell, G.R. (1990), The maximal smoothing principle in density estimation, \emph{Journal of the American Statistical Association}, \bold{85}, 470-477.
}
\author{
T.M. Davies
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
## To be filled
}