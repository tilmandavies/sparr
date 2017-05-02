#' Normal scale (NS) bandwidth selector
#' 
#' Provides the asymptotically optimal isotropic fixed bandwidth for a 2D
#' normal density based on a simple expression.
#' 
#' This function calculates a scalar smoothing bandwidth for kernel density
#' estimates of 2-dimensional data: the optimal value which would minimise the
#' asymptotic mean integrated squared error of the bivariate normal density
#' function, assuming the standard Gaussian kernel function. See Silverman
#' (1986) or Wand and Jones (1995).
#' \describe{
#'   \item{\bold{Effective sample size}}{ The formula requires a sample size,
#'     and this can be minimally
#'     tailored via \code{nstar}. By default, the function simply uses the number
#'     of observations in \code{pp}: \code{nstar = "npoints"}. If \code{pp} is a
#'     \code{\link[spatstat]{ppp.object}} with factor-valued
#'     \code{\link[spatstat]{marks}}, then the user has the option of using
#'     \code{nstar = "geometric"}, which sets the sample size used in the formula
#'     to the geometric mean of the counts of observations of each mark. This can
#'     be useful for e.g. relative risk calculations, see Davies and Hazelton
#'     (2010). Alternatively, the user can specify their own value by simply
#'     supplying a single positive numeric value to \code{nstar}. }
#'   \item{\bold{Spatial scale}}{ For isotropic smoothing, we require a scalar
#'     estimate of the spatial scale of the data. One may use the `robust' estimate
#'     of standard deviation found by a weighted mean of the interquartile ranges
#'     of the \eqn{x}- and \eqn{y}-coordinates of the data respectively
#'     (\code{scaler = "IQR"}). Two other options are the raw mean of the
#'     coordinate-wise standard deviations (\code{scaler = "sd"}), or the square
#'     root of the mean of the two variances (\code{scaler = "var"}). A fourth
#'     option, \code{scaler = "silverman"} (default), sets the scaling constant to
#'     be the minimum of the \code{"IQR"} and \code{"sd"} options; see Silverman
#'     (1986), p. 47. Alternatively, like \code{nstar}, the user can specify their
#'     own value by simply supplying a single positive numeric value to
#'     \code{scaler}.}
#' }
#' 
#' @param pp An object of class \code{\link[spatstat]{ppp}} giving the observed
#'   2D data to be smoothed.
#' @param nstar Optional. Controls the value to use in place of the number of
#'   observations \emph{n} in the normal scale formula. Either a character
#'   string, \code{"npoints"} (default) or \code{"geometric"}, or a positive
#'   numeric value. See `Details'.
#' @param scaler Optional. Controls the value for a scalar representation of
#'   the 2D scale of the data. Either a character string, \code{"silverman"}
#'   (default), \code{"IQR"}, \code{"sd"}, or \code{"var"}; or a positive numeric
#'   value. See `Details'.
#'
#' @return A single numeric value of the estimated bandwidth.
#'
#' @section Warning: The NS bandwidth is an approximation, and assumes
#' \emph{that the target density is bivariate normal}. This is considered rare
#' in most real-world applications. Nevertheless, it remains a quick and easy
#' `rule-of-thumb' method with which one may obtain a smoothing parameter in
#' general applications. Note that a similar expression for the adaptive kernel
#' estimator is not possible (Davies et al., 2017).
#'
#' @author T.M. Davies
#'
#' @references
#' Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel
#' estimation of spatial relative risk, \emph{Statistics in Medicine},
#' \bold{29}(23) 2423-2437.
#'
#' Davies, T.M., Flynn, C.R. and Hazelton, M.L.
#' (2017), On the utility of asymptotic bandwidth selectors for spatially
#' adaptive kernel density estimation, \emph{Submitted}.
#'
#' Silverman, B.W. (1986), \emph{Density Estimation for Statistics and Data Analysis}, Chapman
#' & Hall, New York.
#'
#' Wand, M.P. and Jones, C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.
#'
#' @examples
#' 
#' ## To be filled
#' 
NS <- function(pp, nstar = c("npoints", "geometric"), scaler = c("silverman", "IQR", "sd", "var")){
  if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")
  
  nstar <- processnstar(nstar,pp)
  scaler <- processscaler(scaler,pp)
  
  return(scaler*nstar^(-1/6)) # was: 2*nstar in denominator prior to adjust
}
  
