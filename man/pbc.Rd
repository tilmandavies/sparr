\name{pbc}
\Rdversion{1.1}
\alias{pbc}
\docType{data}
\title{
Primary biliary cirrhosis data
}
\description{
Data of the locations of 761 cases of primary biliary cirrhosis in several adjacent health regions of north-eastern England, along with 3020 controls representing the at-risk population, collected between 1987 and 1994. These data were first presented and analysed by Prince et al. (2001); subsequent analysis of these data in the spirit of \code{\link{sparr}} was performed in Davies and Hazelton (2010). Also included is the polygonal study region.
}
\usage{data(pbc)}
\format{\code{pbc} is a dichotomously marked \code{\link[spatstat]{ppp.object}}, with locations expressed in UK Ordnance Survey Coordinates (km).}

% \describe{
% \tabular{rlll}{
% \code{PBC$x} \tab Easting coordinate (m) \cr
% \code{PBC$y} \tab Northing coordinate (m) \cr
% \code{PBC$marks} \tab Dichotomous \code{\link{factor}}; "case" and "control" \cr
% \code{PBC$window} \tab Geographical study region as a polygonal \code{\link[spatstat]{owin}} with 115 vertices
% }
% }


\source{
Prince et al. (2001), The geographical distribution of primary biliary cirrhosis in a well-defined cohort, \emph{Hepatology}, \bold{34}, 1083-1088.\cr
}
\section{Acknowledgements}{
The authors thank Prof. Peter Diggle at Lancaster University (\url{http://www.lancs.ac.uk/staff/diggle/}) for providing access to these data.
}
\references{
Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.
}
\examples{
data(pbc)
summary(pbc)
plot(pbc)
}

\keyword{datasets}
