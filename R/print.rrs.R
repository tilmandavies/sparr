#' @export
print.rrs <- function(x, ...){
	
	if(all(x$rr>=0)) cat("Relative risk surface\n\n")
	else cat("Log-Relative risk surface\n\n")

	cat("--Numerator (case) density--\n")
	print.bivden(x$f)
	cat("\n--Denominator (control) density--\n")
	print.bivden(x$g)

}