#' @rdname summarysparr
#' @method summary rrs
#' @export
summary.rrs <- function(object, ...){
    if(all(object$rr>=0)) cat("Relative Risk Function.\n\n")
    else cat("Log-Relative Risk Function.\n\n")
    
    if(!all(object$rr>=0)) cat("Estimated risk range [",min(object$rr,na.rm=T),", ",max(object$rr,na.rm=T),"]\n",sep="")
    else cat("Estimated log-risk range\n  [",min(object$rr,na.rm=T),",",max(object$rr,na.rm=T),"]\n",sep="")

    cat("\n--Numerator (case) density--\n")
    summary.bivden(object$f)
    cat("\n--Denominator (control) density--\n")
    summary.bivden(object$g)
}
    
