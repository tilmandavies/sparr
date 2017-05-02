print.bivden <- function(x,...){
	cat("Bivariate kernel density/intensity estimate\n\n")
	if(is.null(x$him)) cat("Fixed bandwidth with h0 =",round(x$h0,4),unitname(x$z)[[2]],"(to 4 d.p.)\n")
	else cat("Adaptive smoothing with h0 =",round(x$h0,4),unitname(x$z)[[2]],"(to 4 d.p.)\n")
	
	cat("No. of observations:",npoints(x$pp),"\n")
}