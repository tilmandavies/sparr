#' Plot tolerance contour classification scheme
#' 
#' Permits illustration of the uniquely identified tolerance contour regions
#' arising from a call to \code{\link{tol.classify}}.
#' 
#' The \code{\link{tol.classify}} function permits identification of 
#' individual significance regions (that is, the tolerance contours). In
#' turn, \code{tol.classplot} may be used to visualise these regions
#' optionally annotated by their unique identification number to better
#' understand the region-specific classifications of the case and control points.
#'  
#' 
#' @param pcpolys A list of polygonal windows, each of class \code{\link[spatstat.geom]{owin}}.
#'  This will almost always be the \code{pcpolys} component of the object
#'  returned by a call to \code{\link{tol.classify}}. 
#' @param add A logical value indicating whether to add the unique regions to an
#'  existing plot (see 'Examples').
#' @param annotate A logical value indicating whether to annotate each unique
#'  region with its identifying number (which will correspond to the uniquely 
#'  split/classified points in a corresponding call to \code{\link{tol.classify}}).
#' @param ... Additional arguments to be passed to \code{\link[graphics]{text}} 
#'  to control the appearance of the annotations when \code{annotate=TRUE}.
#'   
#' @return Plots to the relevant graphics device.
#'
#'
#' @author T. M. Davies
#'
#' @references
#' 
#' Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative
#' risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.
#'
#' Hazelton, M.L. and Davies, T.M. (2009), Inference based on kernel estimates
#' of the relative risk function in geographical epidemiology,
#' \emph{Biometrical Journal}, \bold{51}(1), 98-109.
#'
#' Kelsall, J.E. and Diggle, P.J. (1995), Kernel estimation of relative risk, \emph{Bernoulli},
#' \bold{1}, 3-16.
#'
#' @examples
#' 
#' \dontrun{
#' chrr <- risk(chorley,h0=0.7,tolerate=TRUE)
#' chclass <- tol.classify(chrr,cutoff=0.4)
#' 
#' oldpar <- par(mfrow=c(1,3))
#' #
#' plot(chrr,tol.args=list(levels=0.4))
#' tol.classplot(chclass$pcpolys)
#' plot(Window(chorley))
#' axis(1)
#' axis(2)
#' box(bty="l")
#' tol.classplot(chclass$pcpolys,add=TRUE,col=2,font=2,cex=1.5)
#' #
#' par(oldpar)
#' 
#' }
#' 
#' 
#' 
#' 
#' 
#' @export
tol.classplot <- function(pcpolys,add=FALSE,annotate=TRUE,...){
  if(!is.list(pcpolys)){
    stop("'pcpolys' must be a list")
  }
  
  npol <- length(pcpolys)
  if(npol==0) stop("'pcpolys' is empty")
  
  ocheck <- sapply(pcpolys,is.owin)
  if(!all(ocheck)){
    stop("'pcpolys' must be a list of objects of class 'owin'")
  }
  
  if(!add){
    unwin <- pcpolys[[1]]
    if(npol>1){
      for(i in 2:npol) unwin <- union.owin(unwin,pcpolys[[i]])
    }
    plot(unwin,main="")
  } else {
    for(i in 1:npol) plot(pcpolys[[i]],add=TRUE)
  }
  
  if(annotate){
    cens <- sapply(pcpolys,centroid.owin)
    text(cens[1,],cens[2,],1:npol,...)
  }
}