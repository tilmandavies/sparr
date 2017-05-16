#' Printing sparr objects
#' 
#' \code{print} methods for classes \code{\link{bivden}}, \code{\link{stden}},
#' \code{\link{rrs}}, \code{\link{rrst}} and \code{\link{msden}}.
#' 
#' @aliases print.bivden print.rrs print.msden print.stden print.rrst
#' 
#' @rdname printsparr
#'
#' @param x An object of class \code{\link{bivden}}, \code{\link{stden}},
#'   \code{\link{rrs}}, \code{\link{rrst}}, or \code{\link{msden}}.
#' @param ...  Ignored.
#'
#' @author T.M. Davies
#' @export
print.bivden <- function(x,...){
	cat("Bivariate Kernel Density/Intensity Estimate\n\n")
  if(is.null(x$him)) sm <- "Fixed"
  else sm <- "Adaptive"
  
  cat("Bandwidth\n ",sm,"smoothing with h0 =",round(x$h0,4),unitname(x$z)[[2]],"(to 4 d.p.)\n\n")
  
  cat("No. of observations\n ",npoints(x$pp),"\n")
}
