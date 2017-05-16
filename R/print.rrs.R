#' @rdname printsparr
#' @method print rrs
#' @export
print.rrs <- function(x, ...){
	
	if(all(x$rr>=0)) cat("Relative Risk Surface\n\n")
	else cat("Log-Relative Risk Surface\n\n")

	cat("--Numerator (case) density--\n")
	print.bivden(x$f)
	cat("\n--Denominator (control) density--\n")
	print.bivden(x$g)
}