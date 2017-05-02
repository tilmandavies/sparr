plot.bivden <- function(x, what = c("z", "edge", "bw"), add.pts = FALSE, auto.axes = TRUE, ...){
  
  ellip <- list(...)
  if(is.null(ellip)) ellip <- list()
  
  if(is.null(ellip$main)) ellip$main <- ""
  if(is.null(ellip$box)) ellip$box <- FALSE
  if(is.null(ellip$ribargs)) ellip$ribargs <- list(box=TRUE)
  
  wh <- what[1]
  if(wh=="z"){
    ellip$x <- x$z
    do.call("plot.im",args=ellip)
    if(add.pts) points(x$pp)
  } else if(wh=="edge"){
    ef <- x$q
    if(is.null(ef)){
      ef <- x$z
      ef$v[!is.na(ef$v)] <- 1
      do.call("plot.im",args=ellip)
      if(add.pts) points(x$pp)
      warning("object has no edge correction")
    } else if(!is.im(ef)){
      warning("object has \"diggle\" correction factors as summarised above")
      print(summary(ef))
    } else {
      ellip$x <- ef
      do.call("plot.im",args=ellip)
      if(add.pts) points(x$pp)
    }
  } else if(wh=="bw"){
    ellip$x <- x$him
    if(is.null(ellip$x)){
      ellip$x <- x$z
      ellip$x$v[!is.na(ellip$x$v)] <- x$h0
      warning("object has fixed bandwidth")
    }
    do.call("plot.im",args=ellip)
    if(add.pts) points(x$pp)
  } else {
    stop("invalid 'what'")
  }
  if(auto.axes){
    axis(1)
    axis(2)
    box(bty="l")
    plot(as.polygonal(Window(x$pp)),add=TRUE)
    title(xlab=ellip$xlab,ylab=ellip$ylab)
  }
}