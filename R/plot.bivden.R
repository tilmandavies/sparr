#' Plotting sparr objects
#' 
#' \code{plot} methods for classes \code{\link{bivden}}, \code{\link{stden}},
#' \code{\link{rrs}}, \code{\link{rrst}} and \code{\link{msden}}.
#' 
#' 
#' In all instances, visualisation is deferred to
#' \code{\link[spatstat]{plot.im}}, for which there are a variety of
#' customisations available the user can access via \code{...}. The one
#' exception is when plotting observation-specific \code{"diggle"} edge
#' correction factors---in this instance, a plot of the spatial observations is
#' returned with size proportional to the influence of each correction weight.
#' 
#' When plotting a \code{\link{rrs}} object, a pre-computed \emph{p}-value
#' surface (see argument \code{tolerate} in \code{\link{risk}}) will
#' automatically be superimposed at a significance level of 0.05. Greater
#' flexibility in visualisation is gained by using \code{\link{tolerance}} in
#' conjunction with \code{\link{contour}}.
#' 
#' An \code{\link{msden}}, \code{\link{stden}}, or \code{\link{rrst}} object is plotted as an animation, one pixel image
#' after another, separated by \code{sleep} seconds. If instead you intend the
#' individual images to be plotted in an array of images, you should first set
#' up your plot device layout, and ensure \code{override.par = FALSE} so that
#' the function does not reset these device parameters itself. In such an
#' instance, one might also want to set \code{sleep = 0}.
#' 
#' @aliases plot.bivden plot.rrs plot.msden plot.stden plot.rrst
#'
#' @rdname plotsparr
#'
#'
#' @param x An object of class \code{\link{bivden}}, \code{\link{stden}},
#'   \code{\link{rrs}}, \code{\link{rrst}}, or \code{\link{msden}}.
#' @param what A character string to select plotting of result (\code{"z"};
#'   default); edge-correction surface (\code{"edge"}); or variable bandwidth
#'   surface (\code{"bw"}).
#' @param tselect Either a single numeric value giving the time at which to return the plot, or a vector of length 2 giving an interval of times over which to plot. This argument must respect the stored temporal bound in \code{x$tlim}, else an error will be thrown. By default, the full set of images (i.e. over the entire available time span) is plotted.
#' @param type A character string to select plotting of joint/unconditional spatiotemporal
#'   estimate (default) or conditional spatial density given time.
#' @param fix.range Logical value indicating whether use the same color scale limits for each plot in the sequence. Ignored if the user supplies a pre-defined \code{\link[spatstat]{colourmap}} to the \code{col} argument, which is matched to \code{...} above and passed to \code{\link[spatstat]{plot.im}}. See `Examples'.
#' @param tol.show Logical value indicating whether to show pre-computed tolerance contours on the plot(s). The object \code{x} must already have the relevant \emph{p}-value surface(s) stored in order for this argument to have any effect.
#' @param tol.type A character string used to control the type of tolerance contour displayed; a test for elevated risk (\code{"upper"}), decreased risk (\code{"lower"}), or a two-tailed test (\code{two.sided}).
#' @param tol.args A named list of valid arguments to be passed directly to \code{\link[graphics]{contour}} to control the appearance of plotted contours. Commonly used items are \code{levels}, \code{lty}, \code{lwd} and \code{drawlabels}.
#' @param add.pts Logical value indicating whether to add the observations to
#'   the image plot using default \code{\link{points}}.
#' @param auto.axes Logical value indicating whether to display the plot with
#'   automatically added x-y axes and an `L' box in default styles.
#' @param sleep Single positive numeric value giving the amount of time (in
#'   seconds) to \code{\link[base]{Sys.sleep}} before drawing the next image in
#'   the animation.
#' @param override.par Logical value indicating whether to override the
#'   existing graphics device parameters prior to plotting, resetting
#'   \code{mfrow} and \code{mar}. See `Details' for when you might want to
#'   disable this.
#' @param ...  Additional graphical parameters to be passed to
#'   \code{\link[spatstat]{plot.im}}, or in one instance, to
#'   \code{\link[spatstat]{plot.ppp}} (see `Details').
#'
#' @return Plots to the relevant graphics device.
#'
#' @author T.M. Davies
#'
#' @examples
#' 
#' # to be filled 
#' 
#' @export
plot.bivden <- function(x, what = c("z", "edge", "bw"), add.pts = FALSE, auto.axes = TRUE, sleep = 0.2, override.par = TRUE, ...){
  
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
