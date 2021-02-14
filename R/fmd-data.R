#' Veterinary foot-and-mouth disease outbreak data
#' 
#' Data of the spatial locations and time of farms infected by veterinary foot-and-mouth disease
#' in the county of Cumbria, UK, over a course of nearly 250 days between February and August in 2001.
#' There are 410 infected farms (the cases), and 1866 uninfected farms (the controls). The data
#' have been jittered and randomly thinned by an unspecified amount to preserve anonymity. 
#'  
#' @name fmd
#' @format \code{fmd} is a named list with two members:
#' \describe{
#' \item{\code{$cases}}{
#' An object of class \code{\link[spatstat.geom]{ppp}} giving the spatial locations of the 410 infected
#' farms within a polygonal study region representing the county of Cumbria. The \code{\link[spatstat.geom]{marks}}
#' component of this object contain the integer day of infection (from beginning of study period).
#' }
#' 
#' \item{\code{$controls}}{
#' An object of class \code{\link[spatstat.geom]{ppp}} defined over the same spatial study region with the locations
#' of the 1866 uninfected farms.
#' }
#' }
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @section Acknowledgements: The Animal and Plant Health Agency (APHA), UK, provided permission to use this dataset.
#' 
#' @references 
#' Fernando, W.T.P.S. and Hazelton, M.L. (2014), Generalizing the spatial relative risk function,
#' \emph{Spatial and Spatio-temporal Epidemiology}, \bold{8}, 1-10.
#' 
#' Keeling M, Woolhouse M, Shaw D, Matthews L, Chase-Topping M, Haydon D, et al. (2001),
#' Dynamics of the 2001 UK foot and mouth epidemic: stochastic dispersal in a heterogeneous landscape,
#' \emph{Science}, \bold{294}, 813-817.
#' 
#' Lawson A, Zhou H. (2005), Spatial statistical modeling of disease outbreaks with particular
#' reference to the UK foot and mouth disease (FMD) epidemic of 2001,
#' \emph{Preventative Veterinary Medicine}, \bold{71}, 141-156.
#' 
#' 
#' @examples
#' 
#' data(fmd)
#' summary(fmd$cases)
#' summary(fmd$controls)
#'  
#' par(mfrow=c(1,2))
#' plot(fmd$cases)
#' plot(fmd$controls)
#' 
NULL

#418 farms out of a total of 2813 in the region became infected over the course of the study period,
#with time of infection recorded to the nearest day since the start of the study.
# that remained uninfected over the course of the study period.
# @source \bold{UNSURE IF THESE DATA CAN BE RELEASED WITH sparr}
#\item{\code{$cases.size}}{
#A numeric vector of length 410 giving the size of each infected farm in \code{$cases} as the total animal population.
#}
#\item{\code{$controls.size}}{
#As above, for the uninfected farms.
#}
