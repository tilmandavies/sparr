#' Available global bandwidth range
#' 
#' Gets universally available global bandwidths as represented by several
#' multi-scale density estimate objects
#' 
#' This simple function merely accesses and returns the maximum lower limit and
#' minimum upper limit of all \code{h0range} components of the
#' \code{\link{msden}} objects passed through \code{...}. Natural numeric error
#' arising from any changes to the bandwidth-axis discretisation resolution in
#' the creation of the \code{\link{msden}} objects (i.e. through the
#' `\code{dimz}' argument) means individual global bandwidth ranges can differ
#' slightly between affected multi-scale estimates, even if they are all
#' applied to the same data set. Can additionally be useful when, for example,
#' creating asymmetric relative risk surfaces based on slices of multi-scale
#' densities with respect to the case and control data sets, because the
#' bandwidth factors differ.
#' 
#' Throws an error if one or more of the \code{h0range} components is
#' incompatible (i.e. all \code{h0range} components must overlap).
#' 
#' @param ...  Any number of objects of class \code{\link{msden}}; possibly
#' named.
#'
#' @return A numeric vector of length 2 providing the range of available global
#' bandwidths compatible with all supplied multi-scale density estimates.
#'
#' @author T.M. Davies
#'
#' @seealso \code{\link{multiscale.density}}, \code{\link{multiscale.slice}}
#'
#' @examples
#' 
#' # See ?multiscale.slice
#' 
#' @export
available.h0 <- function(...){
  unpacked <- list(...)
  
  cls <- lapply(unpacked,function(x) inherits(x,"msden"))
  if(!all(unlist(cls))) stop("function arguments must all be of class \"msden\", arising from a call to 'multiscale.density'")
  
  lo <- sapply(unpacked,function(x) x$h0range[1])
  up <- sapply(unpacked,function(x) x$h0range[2])
  rng <- c(max(lo),min(up))
  
  if(rng[1]>=rng[2]) stop("incompatible 'h0range' components -- check bandwidth scales")
  return(rng)
}
