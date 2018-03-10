#' Normal scale (NS) bandwidth selector
#' 
#' Provides the asymptotically optimal fixed bandwidths for spatial or spatiotemporal 
#' normal densities based on a simple expression.
#' 
#' 
#' These functions calculate scalar smoothing bandwidths for kernel density
#' estimates of spatial or spatiotemporal data: the optimal values would minimise the
#' asymptotic mean integrated squared error assuming normally distributed data; see pp. 46-48
#' of Silverman (1986). The \code{NS} function returns a single bandwidth for isotropic smoothing
#' of spatial (2D) data. The \code{NS.spattemp} function returns two values -- one for
#' the spatial margin and another for the temporal margin, based on independently applying
#' the normal scale rule (in 2D and 1D) to the spatial and temporal margins of the supplied data.
#' 
#' \describe{
#'   \item{\bold{Effective sample size}}{ The formula
#'     requires a sample size, and this can be minimally tailored via \code{nstar}.
#'     By default, the function simply uses the number of observations in
#'     \code{pp}: \code{nstar = "npoints"}. Alternatively, the user can specify their own value by simply
#'     supplying a single positive numeric value to \code{nstar}. 
#'     For \code{NS} (not applicable to \code{NS.spattemp}), if \code{pp} is a
#'     \code{\link[spatstat]{ppp.object}} with factor-valued
#'     \code{\link[spatstat]{marks}}, then the user has the option of using
#'     \code{nstar = "geometric"}, which sets the sample size used in the formula
#'     to the geometric mean of the counts of observations of each mark. This can
#'     be useful for e.g. relative risk calculations, see Davies and Hazelton
#'     (2010). 
#'   }
#'   \item{\bold{Spatial (and temporal) scale}}{The \code{scaler} argument is used to specify spatial
#'   (as well as temporal, in use of \code{NS.spattemp}) scale. For isotropic smoothing in the spatial
#'   margin, one may use the `robust' estimate
#'     of standard deviation found by a weighted mean of the interquartile ranges
#'     of the \eqn{x}- and \eqn{y}-coordinates of the data respectively
#'     (\code{scaler = "IQR"}). Two other options are the raw mean of the
#'     coordinate-wise standard deviations (\code{scaler = "sd"}), or the square
#'     root of the mean of the two variances (\code{scaler = "var"}). A fourth
#'     option, \code{scaler = "silverman"} (default), sets the scaling constant to
#'     be the minimum of the \code{"IQR"} and \code{"sd"} options; see Silverman
#'     (1986), p. 47. In use of \code{NS.spattemp} the univariate version of the elected scale
#'     statistic is applied to the recorded times of the data for the temporal bandwidth.
#'     Alternatively, like \code{nstar}, the user can specify their
#'     own value by simply supplying a single positive numeric value to
#'     \code{scaler} for \code{NS}, or a numeric vector of length 2 (in the order of \emph{[<spatial scale>, <temporal scale>]})
#'     for \code{NS.spattemp}.
#'   }
#' }
#' 
#' 
#' @aliases NS.spattemp
#' 
#' @rdname NS
#' 
#' @param pp An object of class \code{\link[spatstat]{ppp}} giving the observed
#'   2D data to be smoothed.
#' @param tt A numeric vector of equal length to the number of points in \code{pp}, 
#' giving the time corresponding to each spatial observation. If unsupplied, 
#' the function attempts to use the values in the \code{\link[spatstat]{marks}} 
#' attribute of the \code{\link[spatstat]{ppp.object}} in \code{pp}.
#' @param nstar Optional. Controls the value to use in place of the number of
#'   observations \emph{n} in the normal scale formula. Either a character
#'   string, \code{"npoints"} (default) or \code{"geometric"} (only possible for \code{NS}), or a positive
#'   numeric value. See `Details'.
#' @param scaler Optional. Controls the value for a scalar representation of
#'   the spatial (and temporal for \code{NS.spattemp}) scale of the data. Either a character string, \code{"silverman"}
#'   (default), \code{"IQR"}, \code{"sd"}, or \code{"var"}; or a positive numeric
#'   value. See `Details'.
#'
#' @return A single numeric value of the estimated spatial bandwidth for \code{NS}, or a named numeric vector of length 2 giving
#' the spatial bandwidth (as \code{h}) and the temporal bandwidth (as \code{lambda}) for \code{NS.spattemp}.
#'
#' @section Warning: The NS bandwidth is an approximation, and assumes
#' \emph{that the target density is normal}. This is considered rare
#' in most real-world applications. Nevertheless, it remains a quick and easy
#' `rule-of-thumb' method with which one may obtain a smoothing parameter. Note that a similar expression for the adaptive kernel
#' estimator is not possible (Davies et al., 2018).
#'
#' @author T.M. Davies
#'
#' @references
#' Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel
#' estimation of spatial relative risk, \emph{Statistics in Medicine},
#' \bold{29}(23) 2423-2437.
#'
#' Davies, T.M., Flynn, C.R. and Hazelton, M.L.
#' (2018), On the utility of asymptotic bandwidth selectors for spatially
#' adaptive kernel density estimation, \emph{Statistics & Probability Letters} [in press].
#'
#' Silverman, B.W. (1986), \emph{Density Estimation for Statistics and Data Analysis}, Chapman
#' & Hall, New York.
#'
#' Wand, M.P. and Jones, C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.
#'
#' @examples
#' 
#' data(pbc)
#' 
#' NS(pbc)
#' NS(pbc,nstar="geometric") # uses case-control marks to replace sample size
#' NS(pbc,scaler="var") # set different scalar measure of spread
#' 
#' data(burk)
#' NS.spattemp(burk$cases)
#' NS.spattemp(burk$cases,scaler="sd") 
#' 
#' @export
NS <- function(pp, nstar = c("npoints", "geometric"), scaler = c("silverman", "IQR", "sd", "var")){
  if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")
  
  nstar <- processnstar(nstar,pp)
  scaler <- processscaler(scaler,pp)
  
  return(scaler*nstar^(-1/6)) # was: 2*nstar in denominator prior to adjust
}
  
