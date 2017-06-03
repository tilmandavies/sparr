#' Burkitt's lymphoma in Uganda
#' 
#' Data of the spatiotemporal locations of Burkitt's lymphoma in the Western Nile district of Uganda
#' from 1960 to 1975.
#' 
#' @name burk
#' @format \code{burk} is a named list with two members:
#' \describe{
#' \item{\code{$cases}}{
#' An object of class \code{\link[spatstat]{ppp}} giving the spatial locations (eastings/northings) 
#' of the 188 cases of Burkitt's lymphoma recorded in individuals of various ages (mostly children); the spatial study region as a polygonal \code{\link[spatstat]{owin}}; as well as the time
#' (in days since 1/1/1960) of each observation stored as the \code{marks} of the points.
#' }
#' 
#' \item{\code{$cases.age}}{
#' A numeric vector of length 188 giving the age of each individual in \code{$cases}.
#' }
#' 
#' }
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @source These data were extracted from the \code{\link[splancs]{burkitt}} object of the \code{splancs} R package;
#' see \cr\cr
#' Rowlingson B. and Diggle P.J. (2017), splancs: Spatial and Space-Time Point Pattern Analysis, R
#' package version 2.01-40; \url{https://CRAN.R-project.org/package=splancs}.
#' 
#' @references 
#' 
#' Bailey T.C. and Gatrell A.C. (1995), \emph{Interactive spatial data analysis}, Longman; Harlow.
#' 
#' @examples
#' data(burk)
#' summary(burk$cases)
#' 
#' par(mfrow=c(1,2))
#' plot(burk$cases)
#' hist(burk$cases.age)
#' 
NULL
