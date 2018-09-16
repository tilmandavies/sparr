#' @rdname plotsparr
#' @method plot rrst
#' @export
plot.rrst <- function(x, tselect = NULL, type = c("joint", "conditional"), fix.range = FALSE, tol.show = TRUE, tol.type = c("upper", "lower", "two.sided"), tol.args = list(levels = 0.05, lty = 1, drawlabels = TRUE), sleep = 0.2, override.par = TRUE, expscale = FALSE, ...){
  ellip <- list(...)
  if(is.null(ellip)) ellip <- list()
  if(is.null(ellip$box)) ellip$box <- FALSE
  if(is.null(ellip$ribargs)) ellip$ribargs <- list(box=TRUE)
  if(!is.null(ellip$zlim)) fix.range <- TRUE
  ellip$log <- FALSE
  mn <- is.null(ellip$main)
  
  
  
  typ <- type[1]
  if(typ=="joint"){
    lst <- x$rr
    plst <- x$P
  } else if(typ=="conditional"){
    lst <- x$rr.cond
    plst <- x$P.cond
  } else stop("Invalid 'type'")
  
  if(override.par) par(mfrow=c(1,1),mar=rep(2,4))
  
  zlimeq <- c(0,min(sapply(lst,max)[sapply(lst,max)>0]))
  zlimconstant <- range(sapply(lst,range))
  
  grt <- as.numeric(names(lst))
  
  if(!is.null(tselect)){
    tsel <- checktsel(tselect)
    if(!all(sapply(tsel,function(y) y>=x$tlim[1]) & sapply(tsel,function(y) y<=x$tlim[2]))) stop(paste("'tselect' must be within allowable time range of",prange(x$tlim)))
    index <- which(grt>=tsel[1]&grt<=tsel[2])
    if(length(index)==0){
      grt <- unique(tsel)
      intrp <- spattemp.slice(x,grt,checkargs=FALSE)
      if(typ=="joint"){
        lst <- intrp$rr
        plst <- intrp$P
      } else {
        lst <- intrp$rr.cond
        plst <- intrp$P.cond
      } 
    } else {
      grt <- grt[index]
      lst <- lst[index]
      plst <- plst[index]
    }
  }
  
  if(!is.null(ellip$zlim)) zlimconstant <- ellip$zlim
  
  if(expscale){
    lst <- lapply(lst,exp)
    ellip$log <- FALSE
    
    if(fix.range&&is.null(ellip$col)){
      # if(is.null(ellip$zlim)){
        ellip$col <- beachcolourmap(range=exp(zlimconstant),sealevel=1)
      # } else {
        # ellip$col <- beachcolourmap(range=ellip$zlim,sealevel=1)
      # }
    }
  }
  
  if(!fix.range) rngs <- lapply(lst,range,na.rm=TRUE)
  
  if(length(lst)==1) sleep <- 0
  
  drawtol <- tol.show&&!is.null(plst)
  if(drawtol){
    plst <- lapply(plst,function(x) t(as.matrix(x)))
    tellip <- tol.args
    tellip$add <- TRUE
    tellip$x <- lst[[1]]$xcol
    tellip$y <- lst[[1]]$yrow
    tol.type <- tol.type[1]
    if(tol.type=="lower"){
      plst <- lapply(plst,function(x) 1-x)
    } else if(tol.type=="two.sided"){
      plst <- lapply(plst,function(x) 2*pmin(x,1-x))
    } else if(tol.type!="upper"){
      stop("invalid 'tol.type'")
    }
  }
  
  for(i in 1:length(lst)){
    dev.hold()
    
    ellip$x <- lst[[i]]
    if(mn) ellip$main <- paste("t =",round(grt[i],5))
    if(diff(range(lst[[i]]))==0&&is.null(ellip$zlim)&&!fix.range) ellip$zlim <- zlimeq
    if(fix.range){
      ellip$zlim <- zlimconstant
      if(expscale) ellip$zlim <- NULL
    } else {
      ellip$zlim <- rngs[[i]]
      if(expscale&&(is.null(ellip$col)||i>1)){
        ellip$col <- beachcolourmap(range=ellip$zlim,sealevel=1)
        ellip$zlim <- NULL
      }
    }
    
    # print(ellip)
    
    do.call("plot.im",ellip)
    
    if(drawtol){
      tellip$z <- plst[[i]]
      suppressWarnings(do.call("contour",tellip))
    }
    
    plot(as.polygonal(Window(x$f$pp)),add=TRUE)
    axis(1)
    axis(2)
    box(bty="l")
    
    dev.flush()
    Sys.sleep(sleep)
  }
  invisible(NULL)
}
