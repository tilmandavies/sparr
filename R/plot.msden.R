#' @rdname plotsparr
#' @method plot msden
#' @export
plot.msden <- function(x, what = c("z", "edge", "bw"), sleep = 0.2, override.par = TRUE, ...){
  wha <- what[1]
  
  ellip <- list(...)
  if(is.null(ellip)) ellip <- list()
  if(is.null(ellip$box)) ellip$box <- FALSE
  if(is.null(ellip$ribargs)) ellip$ribargs <- list(box=TRUE)
  
  if(wha=="z"){
    lst <- x$z
  } else if(wha=="edge"){
    lst <- x$q
    if(is.null(lst)) stop("no edge correction present in multi-scale density object")
  } else if(wha=="bw"){
    lst <- x$him
    if(is.null(ellip$zlim)) ellip$zlim <- range(lapply(lst,range))
  } else {
    stop("invalid 'what'")
  }
  
  if(override.par) par(mfrow=c(1,1),mar=rep(2,4))
  
  hv <- as.numeric(names(lst))
  for(i in 1:length(lst)){
  	dev.hold()
    ellip$x <- lst[[i]]
    ellip$main <- paste("h0 =",round(hv[i],5))
    do.call("plot.im",ellip)
    plot(as.polygonal(Window(x$pp)),add=TRUE)
    axis(1)
    axis(2)
    box(bty="l")
    dev.flush()
    Sys.sleep(sleep)
  }
  invisible(NULL)
}

