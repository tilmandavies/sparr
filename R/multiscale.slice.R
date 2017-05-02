#' Slicing a multi-scale density/intensity object
#' 
#' Takes a slice of the multi-scale density/intensity estimate at a desired
#' global bandwidth
#' 
#' Davies & Baddeley (2017) demonstrate that once a multi-scale
#' density/intensity estimate has been computed, we may take slices parallel to
#' the spatial domain of the trivariate convolution to return the estimate at
#' any desired global bandwidth. This function is the implementation thereof
#' based on a multi-scale estimate resulting from a call to
#' \code{\link{multiscale.density}}.
#' 
#' The function returns an error (unless \code{checkargs = FALSE}) if the
#' requested slice at \code{h0} is not within the available range of
#' pre-computed global bandwidth scalings as defined by the \code{h0range}
#' component of \code{msob}.
#' 
#' Because the contents of the \code{msob} argument, an object of class
#' \code{\link{msden}}, are returned based on a discretised set of global
#' bandwidth scalings, the function internally computes the desired surface as
#' a pixel-by-pixel linear interpolation using the two discretised global
#' bandwidth rescalings that bound the requested \code{h0}. (Thus, numeric
#' accuracy of the slices is improved with an increase to the \code{dimz}
#' argument of the preceding call to \code{multiscale.density} at the cost of
#' additional computing time.)
#' 
#' @param msob An object of class \code{"\link{msden}"} giving the multi-scale
#'   estimate from which to take a slice.
#' @param h0 Desired global bandwidth; the density/intensity estimate
#'   corresponding to which will be returned. This value \bold{must} be in the
#'   available range provided by \code{msob$h0range}; see `Details'.
#' @param checkargs Logical value indicating whether to check validity of
#'   \code{msob} and \code{h0}. Disable only if you know this check will be
#'   unnecessary.
#'
#' @return An object of class \code{\link{bivden}} with components
#' corresponding to the requested slice at \code{h0}. Note that this object
#' will have the component \code{fromms} set to \code{TRUE}.
#'
#' @author T.M. Davies
#'
#' @seealso \code{\link{multiscale.density}}, \code{\link{bivariate.density}}
#'
#' @references
#' Davies, T.M. and Baddeley A. (2017), Fast computation of
#' spatially adaptive kernel estimates, \emph{Submitted}.
#'
#' @export
multiscale.slice <- function(msob,h0,checkargs=TRUE){
  if(checkargs){
    if(!inherits(msob,"msden")) stop("'msob' must be of class \"msden\"")
    h0 <- checkit(h0,"'h0'")
    aran <- msob$h0range
    if(!inside.range(h0,aran)) stop(paste("requested 'h0' outside available range of",prange(aran)))
  }
  
  available.h0 <- as.numeric(names(msob$z))
  zz <- msob$z
  hh <- msob$him
  qq <- msob$q
  
  if(any(available.h0==h0)){
    index <- which(available.h0==h0)
    zres <- zz[[index]]
    hres <- hh[[index]]
    qres <- qq[[index]]
  } else {
    marker <- which(available.h0>h0)[1]
    mindex <- c(marker-1,marker)
    hint <- available.h0[mindex]
    move <- (h0-hint[1])/diff(hint)
    zdiff <- zz[[mindex[2]]]-zz[[mindex[1]]]
    hdiff <- hh[[mindex[2]]]-hh[[mindex[1]]]
    qdiff <- qq[[mindex[2]]]-qq[[mindex[1]]]
    zres <- zz[[mindex[1]]]+move*zdiff
    hres <- hh[[mindex[1]]]+move*hdiff
    if(!is.null(qq)) qres <- qq[[mindex[1]]]+move*qdiff
    else qres <- NULL
  }
  
  result <- list(z=zres,h0=h0,hp=msob$hp,h=msob$h/msob$h0*h0,him=hres,q=qres,gamma=msob$gamma,geometric=msob$geometric,pp=msob$pp,fromms=TRUE)
  class(result) <- "bivden"
  return(result)
}
