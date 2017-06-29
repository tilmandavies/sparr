#' @rdname plotsparr
#' @method plot rrs
#' @export
plot.rrs <- function(x, auto.axes = TRUE, tol.show = TRUE, tol.type = c("upper", "lower", "two.sided"), tol.args = list(levels = 0.05, lty = 1, drawlabels = TRUE), ...){
  ellip <- list(...)
  if(is.null(ellip)) ellip <- list()
  if(is.null(ellip$main)) ellip$main <- ""
  if(is.null(ellip$box)) ellip$box <- FALSE
  if(is.null(ellip$ribargs)) ellip$ribargs <- list(box=TRUE)
  ellip$x <- x$rr
  
  do.call("plot.im",args=ellip)

  if(tol.show&&!is.null(x$P)){
    ps <- t(as.matrix(x$P))
    tellip <- tol.args
    tellip$add <- TRUE
    tellip$x <- x$P$xcol
    tellip$y <- x$P$yrow
    tol.type <- tol.type[1]
    if(tol.type=="lower"){
      ps <- 1-ps
    } else if(tol.type=="two.sided"){
      ps <- 2*pmin(ps,1-ps)
    } else if(tol.type!="upper"){
      stop("invalid 'tol.type'")
    }
    tellip$z <- ps
    suppressWarnings(do.call("contour",tellip))
  }
  
  plot(as.polygonal(Window(x$f$pp)),add=TRUE)
  if(auto.axes){
    axis(1)
    axis(2)
    box(bty="l")
    title(xlab=ellip$xlab,ylab=ellip$ylab)
  }
}