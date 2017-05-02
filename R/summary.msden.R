#' Summarising an estimated multi-scale density object
#' 
#' \code{print} and \code{summary} methods for class \code{\link{msden}}
#' 
#' 
#' @aliases summary.msden print.msden
#' @param x,object An object of class \code{\link{msden}} resulting from a call
#' to \code{\link{multiscale.density}}.
#' @param ...  Ignored.
#' @author T.M. Davies
summary.msden <- function(object,...){
  cat("Multi-scale bivariate adaptive kernel density/intensity estimate\n\n")
  h0r <- round(object$h0range,4)
  h0v <- as.numeric(names(object$z))
  h0l <- length(h0v)
  cat("Available global bandwidth range (",h0r[1],", ",h0r[2],") ",unitname(object$z[[1]])[[2]],";\n",sep="")
  cat("  discretised sequence of length",h0l,"\b: ")
  cat(round(h0v,4),"\n",sep=", ")
  cat("\b\b\b.")
  if(!is.null(object$q)) cat("\n\nEdge-corrected.\n")
  cat("\nNo. of observations:",npoints(object$pp),"\n")
}
