#' Burkitt's lymphoma in Uganda
#' 
#' Data of the spatiotemporal locations of Burkitt's lymphoma in the Western Nile district of Uganda
#' from 1960 to 1975.
#' 
#' @name burk
#' @format \code{burk} is a named list with three members:
#' \describe{
#' \item{\code{$cases}}{
#' An object of class \code{\link[spatstat.geom]{ppp}} giving the spatial locations (eastings/northings) 
#' of the 188 cases of Burkitt's lymphoma recorded in individuals of various ages (mostly children); the spatial study region as a polygonal \code{\link[spatstat.geom]{owin}}; as well as the time
#' (in days since 1/1/1960) of each observation stored as the \code{marks} of the points.
#' }
#' 
#' \item{\code{$cases.age}}{
#' A numeric vector of length 188 giving the age of each individual in \code{$cases}.
#' }
#' 
#' \item{\code{$controls}}{
#' An object of class \code{\link[spatstat.geom]{ppp}} giving 500 \bold{artificially simulated} spatial-only
#' observations to pose as a `control' data set representing the at-risk population. The data were
#' generated from a smooth kernel estimate of the spatial margin of the cases. The similarity between the case point distribution
#' and the true at-risk population dispersion can be seen in e.g. Figure 2 of Middleton and Greenland (1954).
#' 
#' }
#' 
#' }
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @source The case data were extracted from the \code{\link[splancs]{burkitt}} object of the \code{splancs} R package;
#' see \cr\cr
#' Rowlingson B. and Diggle P.J. (2017), splancs: Spatial and Space-Time Point Pattern Analysis, R
#' package version 2.01-40; \url{https://CRAN.R-project.org/package=splancs}.
#' 
#' @references 
#' 
#' Bailey, T.C. and Gatrell, A.C. (1995), \emph{Interactive spatial data analysis}, Longman; Harlow.
#' 
#' Middleton, J.F.M. and Greenland, D.J. (1954), Land and population in West Nile District, Uganda, \emph{The Geographical Journal}, \bold{120}, 446--455.
#' 
#' @examples
#' data(burk)
#' summary(burk$cases)
#' 
#' par(mfrow=c(1,3))
#' plot(burk$cases)
#' plot(burk$controls)
#' plot(density(marks(burk$cases)),xlim=range(marks(burk$cases)))
NULL
