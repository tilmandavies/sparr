#' @rdname printsparr
#' @method print msden
#' @export
print.msden <- function(x,...){
  cat("Multi-scale Adaptive Kernel Density/Intensity Estimate\n\n")
  cat("Available global bandwidth range\n  (",round(x$h0range[1],4),", ",round(x$h0range[2],4),") ",unitname(x$z[[1]])[[2]]," (to 4 d.p.)\n\n",sep="")
  cat("No. of observations\n ",npoints(x$pp),"\n")
}