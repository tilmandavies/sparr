print.msden <- function(x,...){
  cat("Multi-scale bivariate adaptive kernel density/intensity estimate\n\n")
  cat("Available global bandwidth range (",round(x$h0range[1],4),", ",round(x$h0range[2],4),") ",unitname(x$z[[1]])[[2]]," (to 4 d.p.)\n",sep="")
  cat("No. of observations:",npoints(x$pp),"\n")
}