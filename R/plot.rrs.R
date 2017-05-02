#' @export
plot.rrs <- function(x, auto.axes = TRUE, ...){
  ellip <- list(...)
  if(is.null(ellip$main)) ellip$main <- ""
  if(is.null(ellip$box)) ellip$box <- FALSE
  if(is.null(ellip$ribargs)) ellip$ribargs <- list(box=TRUE)
  ellip$x <- x$rr
  
  do.call("plot.im",args=ellip)
  
  if(!is.null(x$P)) contour(x$rr$xcol,x$rr$yrow,t(as.matrix(x$P)),level=0.05,drawlabels=FALSE,add=TRUE)
  
  if(auto.axes){
    axis(1)
    axis(2)
    box(bty="l")
    plot(Window(x$f$pp),add=TRUE)
    title(xlab=ellip$xlab,ylab=ellip$ylab)
  }
}