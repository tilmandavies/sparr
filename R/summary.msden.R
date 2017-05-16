#' @rdname summarysparr
#' @method summary msden
#' @export
summary.msden <- function(object,...){
  cat("Multi-scale Adaptive Kernel Density/Intensity Estimate\n\n")
  h0r <- round(object$h0range,4)
  h0v <- as.numeric(names(object$z))
  h0l <- length(h0v)
  cat("Available global bandwidth range\n  (",h0r[1],", ",h0r[2],") ",unitname(object$z[[1]])[[2]],"\n",sep="")
  cat("  Discretised sequence of length",h0l,"\b: ")
  cat(round(h0v,4),"\n",sep=", ")
  cat("\b\b\b.")
  cat("\n\nNo. of observations\n ",npoints(object$pp),"\n")
  
  cat("\nEvaluation per slice\n  ",nrow(object$z[[1]])," x ",ncol(object$z[[1]])," rectangular grids\n  ",sum(!is.na(as.vector(as.matrix(object$z[[1]]))))," grid cells out of ",prod(dim(object$z[[1]]))," fall inside study region\n",sep="")
}
