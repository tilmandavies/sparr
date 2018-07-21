#' @rdname plotsparr
#' @method plot stden
#' @export
plot.stden <- function(x, tselect = NULL, type = c("joint", "conditional"), fix.range = FALSE, sleep = 0.2, override.par = TRUE, ...){
  ellip <- list(...)
  if(is.null(ellip)) ellip <- list()
  if(is.null(ellip$box)) ellip$box <- FALSE
  if(is.null(ellip$ribargs)) ellip$ribargs <- list(box=TRUE)
  if(is.null(ellip$log)) ellip$log <- FALSE
  mn <- is.null(ellip$main)

  typ <- type[1]
  if(typ=="joint"){
    lst <- x$z
  } else if(typ=="conditional"){
    lst <- x$z.cond
  } else stop("Invalid 'type'")
  
  if(override.par) par(mfrow=c(1,1),mar=rep(2,4))
  
  zlimeq <- c(0,min(sapply(lst,max)[sapply(lst,max)>0]))
  zlimconstant <- range(sapply(lst,range))
  if(ellip$log&&fix.range) zlimconstant <- log(zlimconstant)
  grt <- as.numeric(names(lst))
  
  if(!is.null(tselect)){
    tsel <- checktsel(tselect)
    if(!all(sapply(tsel,function(y) y>=x$tlim[1]) & sapply(tsel,function(y) y<=x$tlim[2]))) stop(paste("'tselect' must be within allowable time range of",prange(x$tlim)))
    index <- which(grt>=tsel[1]&grt<=tsel[2])
    if(length(index)==0){
      grt <- unique(tsel)
      intrp <- spattemp.slice(x,grt,checkargs=FALSE)
      if(typ=="joint") lst <- intrp$z
      else lst <- intrp$z.cond
    } else {
      grt <- grt[index]
      lst <- lst[index]
    }
  }
  
  if(length(lst)==1) sleep <- 0
  
  for(i in 1:length(lst)){
    dev.hold()
    ellip$x <- lst[[i]]
    if(mn) ellip$main <- paste("t =",round(grt[i],5))
    if(diff(range(lst[[i]]))==0&&is.null(ellip$zlim)&&!fix.range) ellip$zlim <- zlimeq
    if(fix.range) ellip$zlim <- zlimconstant
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

