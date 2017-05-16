#' Summarising sparr objects
#' 
#' \code{summary} methods for classes \code{\link{bivden}}, \code{\link{stden}},
#' \code{\link{rrs}}, \code{\link{rrst}} and \code{\link{msden}}.
#' 
#' @aliases summary.bivden summary.rrs summary.msden summary.stden summary.rrst
#' 
#' @rdname summarysparr
#'
#' @param object An object of class \code{\link{bivden}}, \code{\link{stden}},
#'   \code{\link{rrs}}, \code{\link{rrst}}, or \code{\link{msden}}.
#' @param ...  Ignored.
#'
#' @author T.M. Davies
#' @export
summary.bivden <- function(object, ...){

	print.bivden(x=object)
	
  W <- Window(object$pp)
  wt <- summary(W)$type
  wx <- W$xrange
  wy <- W$yrange
  cat("\nSpatial bound\n  Type: ",wt,"\n  2D enclosure: [",wx[1],", ",wx[2],"] x [",wy[1],", ",wy[2],"]\n",sep="")
  
  cat("\nEvaluation\n  ",nrow(object$z)," x ",ncol(object$z)," rectangular grid\n  ",sum(!is.na(as.vector(as.matrix(object$z))))," grid cells out of ",prod(dim(object$z))," fall inside study region\n",sep="")
	
	cat("  Density/intensity range [",min(object$z,na.rm=TRUE),", ",max(object$z,na.rm=TRUE),"]\n",sep="")
}
