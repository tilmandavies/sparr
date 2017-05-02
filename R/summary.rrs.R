#' Summarising an estimated relative risk function object
#' 
#' \code{print} and \code{summary} methods for class \code{\link{rrs}}.
#' 
#' 
#' @aliases summary.rrs print.rrs
#' @param x,object An object of class \code{\link{rrs}} resulting from a call
#' to \code{\link{risk}}.
#' @param ...  Ignored.
#' @author T.M. Davies
summary.rrs <- function(object, ...){
    if(all(object$rr>=0)) cat("Relative risk function.\n\n")
    else cat("Log-Relative risk function.\n\n")
    
    
    if(!all(object$rr>=0)) cat("Estimated risk range [",min(object$rr,na.rm=T),",",max(object$rr,na.rm=T),"]\n",sep="")
    else cat("Estimated log-risk range [",min(object$rr,na.rm=T),",",max(object$rr,na.rm=T),"]\n",sep="")
#    cat(sum(!is.na(as.vector(object$rsM))),"grid cells out of",prod(dim(object$rsM)),"fall inside study region.\n")
    
    #cat("Surface (Z) summary:\n")
#    print(summary(as.vector(object$rsM)))
#    
    cat("\n--Numerator (case) density--\n")
    summary.bivden(object$f)
    cat("\n--Denominator (control) density--\n")
    summary.bivden(object$g)
    
}
    
