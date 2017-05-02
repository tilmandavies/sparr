\name{available.h0}
\alias{available.h0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Available global bandwidth range
}
\description{
Gets universally available global bandwidths as represented by several multi-scale density estimate objects
}
\usage{
available.h0(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{
Any number of objects of class \code{\link{msden}}; possibly named.
}
}
\details{
This simple function merely accesses and returns the maximum lower limit and minimum upper limit of all \code{h0range} components of the \code{\link{msden}} objects passed through \code{...}. Natural numeric error arising from any changes to the bandwidth-axis discretisation resolution in the creation of the \code{\link{msden}} objects (i.e. through the `\code{dimz}' argument) means individual global bandwidth ranges can differ slightly between affected multi-scale estimates, even if they are all applied to the same data set. Can additionally be useful when, for example, creating asymmetric relative risk surfaces based on slices of multi-scale densities with respect to the case and control data sets, because the bandwidth factors differ.

Throws an error if one or more of the \code{h0range} components is incompatible (i.e. all \code{h0range} components must overlap).
}
\value{
A numeric vector of length 2 providing the range of available global bandwidths compatible with all supplied multi-scale density estimates.
}
\author{
T.M. Davies
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{multiscale.density}}, \code{\link{multiscale.slice}}
}

\examples{
# To be filled.
}
