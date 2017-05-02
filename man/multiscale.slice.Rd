\name{multiscale.slice}
\alias{multiscale.slice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Slicing a multi-scale density/intensity object
}
\description{
Takes a slice of the multi-scale density/intensity estimate at a desired global bandwidth
}
\usage{
multiscale.slice(msob, h0, checkargs = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{msob}{
An object of class \code{"\link{msden}"} giving the multi-scale estimate from which to take a slice.
}
  \item{h0}{
Desired global bandwidth; the density/intensity estimate corresponding to which will be returned. This value \bold{must} be in the available range provided by \code{msob$h0range}; see `Details'.
}
  \item{checkargs}{
Logical value indicating whether to check validity of \code{msob} and \code{h0}. Disable only if you know this check will be unnecessary.
}
}
\details{
Davies & Baddeley (2017) demonstrate that once a multi-scale density/intensity estimate has been computed, we may take slices parallel to the spatial domain of the trivariate convolution to return the estimate at any desired global bandwidth. This function is the implementation thereof based on a multi-scale estimate resulting from a call to \code{\link{multiscale.density}}.

The function returns an error (unless \code{checkargs = FALSE}) if the requested slice at \code{h0} is not within the available range of pre-computed global bandwidth scalings as defined by the \code{h0range} component of \code{msob}.

Because the contents of the \code{msob} argument, an object of class \code{\link{msden}}, are returned based on a discretised set of global bandwidth scalings, the function internally computes the desired surface as a pixel-by-pixel linear interpolation using the two discretised global bandwidth rescalings that bound the requested \code{h0}. (Thus, numeric accuracy of the slices is improved with an increase to the \code{dimz} argument of the preceding call to \code{multiscale.density} at the cost of additional computing time.)

}
\value{
An object of class \code{\link{bivden}} with components corresponding to the requested slice at \code{h0}. Note that this object will have the component \code{fromms} set to \code{TRUE}.
}
\references{
Davies, T.M. and Baddeley A. (2017), Fast computation of spatially adaptive kernel estimates, \emph{Submitted}.\cr\cr
}
\author{
T.M. Davies
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{multiscale.density}}, \code{\link{bivariate.density}}
}

%\examples{ }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
