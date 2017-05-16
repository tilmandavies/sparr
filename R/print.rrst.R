#' @rdname printsparr
#' @method print rrst
#' @export
print.rrst <- function(x, ...){

  if(all(sapply(x$rrc,min,na.rm=TRUE)>=0)) cat("Spatiotemporal Relative Risk Surface\n\n")
  else cat("Spatiotemporal Log-Relative Risk Surface\n\n")
  
  cat("--Numerator (case) density--\n")
  print(x$f)
  cat("\n--Denominator (control) density--\n")
  print(x$g)
}