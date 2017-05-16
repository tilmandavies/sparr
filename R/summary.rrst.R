#' @rdname summarysparr
#' @method summary rrst
#' @export
summary.rrst <- function(object, ...){
  
  if(all(sapply(object$rrc,min,na.rm=TRUE)>=0)) cat("Spatiotemporal Relative Risk Surface\n\n")
  else cat("Spatiotemporal Log-Relative Risk Surface\n\n")
  
  cat("--Numerator (case) density--\n")
  summary(object$f)
  cat("\n--Denominator (control) density--\n")
  summary(object$g)
}