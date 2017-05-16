#' @rdname printsparr
#' @method print stden
#' @export
print.stden <- function(x,...){
  cat("Spatiotemporal Kernel Density Estimate\n\n")
  cat("Bandwidths\n  h =",round(x$h,4),"(spatial)\n  lambda =",round(x$lambda,4),"(temporal)\n\n")
  cat("No. of observations\n ",npoints(x$pp),"\n")
}