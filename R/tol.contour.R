#' Plot tolerance contours
#' 
#' Draw contours based on a \emph{p}-value matrix.
#' 
#' Note that no checks on the numeric content of \code{pim} are made. The
#' function assumes the pixel \code{\link[spatstat]{im}}age of \emph{p}-values
#' in \code{pim} is supplied with respect to an upper-tailed test for elevated
#' risk (this is exactly the way the \emph{p}-value surface is returned when
#' \code{\link{tolerance}} is used). This is important if one makes subsequent
#' use of \code{test}, which manipulates the \emph{p}-values to draw at desired
#' significance \code{levels}.
#' 
#' @param pim A pixel \code{\link[spatstat]{im}}age of \emph{p}-values,
#'   typically obtained from a call to \code{\link{tolerance}}, computed with
#'   respect to a test for elevated risk.
#' @param test An optional character string giving the type of manipulation to
#'   be applied to the \emph{p}-values, corresponding to a test for significantly
#'   elevated risk (\code{"upper"}; default); for reduced risk (\code{"lower"});
#'   or for both (\code{"two-sided"}).
#' @param ...  Additional arguments to be passed to \code{\link{contour}}.
#'   Commonly used options include \code{add} (to superimpose the contours upon
#'   an existing plot); \code{levels} (to control the specific significance
#'   levels at which to delineate the \emph{p}-values); and \code{lty} or
#'   \code{lwd} for aesthetics.
#'
#' @return Opens a new graphics device and displays a \code{\link{contour}}
#' plot if \code{add = FALSE}, otherwise adds the contours to the plot in the
#' existing active graphics device.
#'
#' @author T. M. Davies
#'
#' @examples
#' 
#' # See ?tolerance
#' 
#' @export
tol.contour <- function(pim, test = c("upper", "lower", "two-sided"), ...){
  if(!inherits(pim,"im")){
    stop("'pim' must be an object of class 'im', typically arising from a call to 'tolerance'")
  }
  tt <- t(as.matrix(pim))
  test <- test[1]
  if(test=="lower"){
    tt <- 1-tt
  } else if(test=="two-sided"){
    tt <- 2*pmin(tt,1-tt)
  } else if(test!="upper"){
    stop("invalid 'test'")
  }
  contour(x=pim$xcol,y=pim$yrow,z=tt,...)
}
