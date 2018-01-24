#' Random point generation inside polygon
#' 
#' Generates a random point pattern of \eqn{n} iid points with any specified distribution based on a pixel image and a corresponding polygonal window.
#' 
#' This function is a deliberate variant of \code{\link[spatstat]{rpoint}} (Baddeley et. al, 2015), to be accessed when the user desires 
#' a randomly generated point pattern based on a pixel image, but wants the window of the point pattern to be a corresponding irregular polygon, as opposed to a binary
#' image mask (which, when converted to a polygon directly, gives jagged edges based on the union of the pixels). When the user specifies their own polygonal window, a \code{while} loop is called and repeated as many
#' times as necessary (up to \code{maxpass} times) to find \code{n} points inside \code{w} (when \code{w = NULL}, then the aforementioned union of the pixels of \code{z}
#' is used, obtained via \code{as.polygonal(Window(z))}). The loop is necessary because the standard behaviour of \code{\link[spatstat]{rpoint}} can (and often does)
#' yield points that sit in corners of pixels which lie outside the corresponding \code{w}.
#' 
#' The \code{correction} argument is used to determine how many points are generated initially,
#' which will be \code{ceiling(correction*n)}; to minimise the number of required passes over the loop this is by default set to give a number slightly higher than the requested \code{n}.
#' 
#' An error is thrown if \code{Window(z)} and \code{w} do not overlap.
#' 
#' @aliases rimpoly
#' 
#' @rdname rimpoly
#' 
#' @param n Number of points to generate.
#' @param z A pixel image of class \code{\link[spatstat]{im}} defining the probability density of the points, possibly unnormalised.
#' @param w A polygonal window of class \code{\link[spatstat]{owin}}. See `Details'.
#' @param correction An adjustment to the number of points generated at the initial pass of the internal loop in an effort to minimise the total number of passes required to reach \eqn{n} points. See `Details'.
#' @param maxpass The maximum number of passes allowed before the function exits. If this is reached before \eqn{n} points are found that fall within \code{w}, a warning is issued.
#'
#' @return An object of class \code{\link[spatstat]{ppp}} containing the \code{n} generated points, defined with the polygonal \code{\link[spatstat]{owin}}, \code{w}.
#'
#' @author T.M. Davies
#'
#' @references
#' Baddeley, A., Rubak, E. and Turner, R. (2015) \emph{Spatial Point Patterns: Methodology and Applications with R}, Chapman and Hall/CRC Press, UK.
#'
#' @examples
#' 
#' data(pbc)
#' Y <- bivariate.density(pbc,h0=2.5,res=25)
#' 
#' # Direct use of 'rpoint':
#' A <- rpoint(500,Y$z)
#' npoints(A)
#' 
#' # Using 'rimpoly' without supplying polygon:
#' B <- rimpoly(500,Y$z)
#' npoints(B)
#' 
#' # Using 'rimpoly' with the original pbc polygonal window:
#' C <- rimpoly(500,Y$z,Window(Y$pp))
#' npoints(C)
#' 
#' par(mfrow=c(1,3))
#' plot(A,main="rpoint")
#' plot(B,main="rimpoly (no polygon supplied)")
#' plot(C,main="rimpoly (original polygon supplied)")
#' 
#' @export
rimpoly <- function(n,z,w=NULL,correction=1.1,maxpass=50){
  if(!is.im(z)) stop("'z' must be a pixel image (spatstat class \"im\")")
  if(is.null(w)) w <- as.polygonal(Window(z))
  if(is.null(intersect.owin(Window(z),w,fatal=FALSE))) stop("'z' and 'w' must overlap")
  
  genblock <- ceiling(correction*n)
  pass <- 1
  result <- matrix(NA,1,2)
  while(((nrow(result)-1)<n)&&(pass<maxpass)){
    tempr <- rpoint(ceiling(genblock/pass),z)
    tempp <- suppressWarnings(ppp(tempr$x,tempr$y,window=w))
    result <- rbind(result,cbind(tempp$x,tempp$y))
    pass <- pass+1
  }
  
  result <- result[2:(n+1),]
  if((pass==maxpass)&&(nrow(result)<n)) warning(paste("Maximum number of passes (",maxpass,") reached with only ",nrow(result)," points generated inside polygon",sep=""))
  
  return(ppp(x=result[,1],y=result[,2],window=w))
}