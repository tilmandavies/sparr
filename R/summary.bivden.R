#' Summarising an estimated bivariate density/intensity object
#' 
#' \code{print} and \code{summary} methods for class \code{\link{bivden}}.
#' 
#' 
#' @aliases summary.bivden print.bivden
#' @param x,object An object of class \code{\link{bivden}} resulting from a
#' call to \code{\link{bivariate.density}}.
#' @param ...  Ignored.
#' @author T.M. Davies
summary.bivden <- function(object, ...){

	print.bivden(x=object)
	
	cat("\nEvaluated over a",nrow(object$z),"by",ncol(object$z),"rectangular grid.\n",sum(!is.na(as.vector(as.matrix(object$z)))),"grid cells out of",prod(dim(object$z)),"fall inside study region.\n\n")
	
	#cat("Estimated density description\n")
	
	cat("Estimated density/intensity range [",min(object$z,na.rm=TRUE),",",max(object$z,na.rm=TRUE),"].\n",sep="")
	
	#print(summary(as.vector(object$Zm)))

}
	
