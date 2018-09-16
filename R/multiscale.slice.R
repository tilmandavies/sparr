#' Slicing a multi-scale density/intensity object
#' 
#' Takes slices of a multi-scale density/intensity estimate at desired
#' global bandwidths
#' 
#' Davies & Baddeley (2018) demonstrate that once a multi-scale
#' density/intensity estimate has been computed, we may take slices parallel to
#' the spatial domain of the trivariate convolution to return the estimate at
#' any desired global bandwidth. This function is the implementation thereof
#' based on a multi-scale estimate resulting from a call to
#' \code{\link{multiscale.density}}.
#' 
#' The function returns an error if the
#' requested slices at \code{h0} are not all within the available range of
#' pre-computed global bandwidth scalings as defined by the \code{h0range}
#' component of \code{msob}.
#' 
#' Because the contents of the \code{msob} argument, an object of class
#' \code{\link{msden}}, are returned based on a discretised set of global
#' bandwidth scalings, the function internally computes the desired surface as
#' a pixel-by-pixel linear interpolation using the two discretised global
#' bandwidth rescalings that bound each requested \code{h0}. (Thus, numeric
#' accuracy of the slices is improved with an increase to the \code{dimz}
#' argument of the preceding call to \code{multiscale.density} at the cost of
#' additional computing time.)
#' 
#' @param msob An object of class \code{\link{msden}} giving the multi-scale
#'   estimate from which to take slices.
#' @param h0 Desired global bandwidth(s); the density/intensity estimate
#'   corresponding to which will be returned. A numeric vector. All values \bold{must} be in the
#'   available range provided by \code{msob$h0range}; see `Details'.
#' @param checkargs Logical value indicating whether to check validity of
#'   \code{msob} and \code{h0}. Disable only if you know this check will be
#'   unnecessary.
#'
#' @return If \code{h0} is scalar, an object of class \code{\link{bivden}} with components
#' corresponding to the requested slice at \code{h0}. If \code{h0} is a vector, a list of objects
#' of class \code{\link{bivden}}. 
#' 
#' @author T.M. Davies
#'
#' @seealso \code{\link{multiscale.density}}, \code{\link{bivariate.density}}
#'
#' @references
#' Davies, T.M. and Baddeley A. (2018), Fast computation of
#' spatially adaptive kernel estimates, \emph{Statistics and Computing}, \bold{28}(4), 937-956.
#' 
#' @examples
#' \donttest{
#' data(chorley) # Chorley-Ribble data (package 'spatstat')
#' ch.multi <- multiscale.density(chorley,h0=1,h0fac=c(0.5,2))
#' 
#' available.h0(ch.multi)
#' ch.slices <- multiscale.slice(ch.multi,h0=c(0.7,1.1,1.6))
#' 
#' par(mfcol=c(2,3)) # plot each density and edge-correction surface
#' for(i in 1:3) { plot(ch.slices[[i]]$z); plot(ch.slices[[i]]$q) }
#' }
#' 
#' @export
multiscale.slice <- function(msob,h0,checkargs=TRUE){
  if(checkargs){
    if(!inherits(msob,"msden")) stop("'msob' must be of class \"msden\"")
    if(!is.vector(h0)||!is.numeric(h0)) stop("'h0' must be a numeric vector")
    
    h0 <- sapply(h0,checkit,str="'h0'")
    aran <- msob$h0range
    
    if(!all(sapply(h0,function(x) x>=aran[1]) & sapply(h0,function(x) x<=aran[2]))) stop(paste("at least one requested 'h0' is outside available range of",prange(aran)))
    # was:
    # if(!inside.range(h0,aran)) stop(paste("requested 'h0' outside available range of",prange(aran)))
  }
  
  avail <- names(msob$z)
  zz <- msob$z
  hh <- msob$him
  qq <- msob$q

  hlen <- length(h0)
  if(hlen==1){
    slc <- ms.slice.single(h0,avail,zz,hh,qq)
    result <- list(z=slc$z,h0=h0,hp=msob$hp,h=msob$h/msob$h0*h0,him=slc$h,q=slc$q,gamma=msob$gamma,geometric=msob$geometric,pp=msob$pp) 
    class(result) <- "bivden"
  } else {
    result <- list()
    for(i in 1:hlen){
      slc <- ms.slice.single(h0[i],avail,zz,hh,qq)
      retob <- list(z=slc$z,h0=h0[i],hp=msob$hp,h=msob$h/msob$h0[i]*h0[i],him=slc$h,q=slc$q,gamma=msob$gamma,geometric=msob$geometric,pp=msob$pp) 
      class(retob) <- "bivden"
      result[[i]] <- retob
    }
  }
  
  return(result)
}

ms.slice.single <- function(V,avail,zz,hh,qq,warn=FALSE){
  la <- length(avail)
  if(any(avail==as.character(V))){
    index <- which(avail==as.character(V))
    zres <- zz[[index]]
    hres <- hh[[index]]
    qres <- qq[[index]]
  } else {
    marker <- as.numeric(avail)>V
    if(sum(marker)==la){
      zres <- zz[[1]]
      hres <- hh[[1]]
      qres <- qq[[1]]
      if(warn) warning("lower index mismatch")
    } else if(sum(marker)==0){
      zres <- zz[[la]]
      hres <- hh[[la]]
      qres <- qq[[la]]
      if(warn) warning("upper index mismatch")
    } else {
      marker <- which(marker)[1]
      mindex <- c(marker-1,marker)
      hint <- as.numeric(avail)[mindex]
      move <- (V-hint[1])/diff(hint)
      zdiff <- zz[[mindex[2]]]-zz[[mindex[1]]]
      hdiff <- hh[[mindex[2]]]-hh[[mindex[1]]]
      qdiff <- qq[[mindex[2]]]-qq[[mindex[1]]]
      zres <- zz[[mindex[1]]]+move*zdiff
      hres <- hh[[mindex[1]]]+move*hdiff
      if(!is.null(qq)) qres <- qq[[mindex[1]]]+move*qdiff
      else qres <- NULL
    }
  }
  return(list(z=zres,h=hres,q=qres))
}
