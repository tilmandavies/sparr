#' Printing a bivariate density/intensity, relative risk surface, or
#' multi-scale density/intensity
#' 
#' \code{print} methods for classes \code{"\link{bivden}"}, \code{"\link{rrs}"},
#' and \code{"\link{msden}"}.
#' 
#' @aliases print.rrs print.msden
#'
#' @param x An object of class \code{\link{bivden}}, \code{\link{rrs}}, or
#'   \code{\link{msden}}.
#' @param ...  Ignored.
#'
#' @author T.M. Davies
#' @export
print.bivden <- function(x,...){
	cat("Bivariate kernel density/intensity estimate\n\n")
	if(is.null(x$him)) cat("Fixed bandwidth with h0 =",round(x$h0,4),unitname(x$z)[[2]],"(to 4 d.p.)\n")
	else cat("Adaptive smoothing with h0 =",round(x$h0,4),unitname(x$z)[[2]],"(to 4 d.p.)\n")
	
	cat("No. of observations:",npoints(x$pp),"\n")
}