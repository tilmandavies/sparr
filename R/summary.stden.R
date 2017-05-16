#' @rdname summarysparr
#' @method summary stden
#' @export
summary.stden <- function(object, ...){
  
  print.stden(x=object)
  
  W <- Window(object$pp)
  wt <- summary(W)$type
  wx <- W$xrange
  wy <- W$yrange
  cat("\nSpatial bound\n  Type: ",wt,"\n  2D enclosure: [",wx[1],", ",wx[2],"] x [",wy[1],", ",wy[2],"]\n",sep="")
  
  cat("\nTemporal bound\n  [",object$tlim[1],", ",object$tlim[2],"]\n",sep="")
  
  cat("\nEvaluation\n  ",nrow(object$z[[1]])," x ",ncol(object$z[[1]])," x ",length(object$z)," trivariate lattice\n",sep="")
  
  minden <- min(sapply(object$z,min,na.rm=TRUE))
  maxden <- max(sapply(object$z,max,na.rm=TRUE))
  cat("  Density range: [",minden,", ",maxden,"]\n",sep="")
}
