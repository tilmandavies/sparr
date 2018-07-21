#' Slicing a spatiotemporal object
#' 
#' Takes slices of the spatiotemporal kernel density or 
#' relative risk function estimate at desired times
#' 
#' 
#' Contents of the \code{stob} argument are returned based on a discretised set of times.
#' This function internally computes the desired surfaces as
#' pixel-by-pixel linear interpolations using the two discretised times
#' that bound each requested \code{tt}.
#' 
#' The function returns an error if any of the
#' requested slices at \code{tt} are not within the available range of
#' times as given by the \code{tlim}
#' component of \code{stob}.
#' 
#' @param stob An object of class \code{\link{stden}} or \code{\link{rrst}} giving the spatiotemporal
#'   estimate from which to take slices.
#' @param tt Desired time(s); the density/risk surface estimate
#'   corresponding to which will be returned. This value \bold{must} be in the
#'   available range provided by \code{stob$tlim}; see `Details'.
#' @param checkargs Logical value indicating whether to check validity of
#'   \code{stob} and \code{tt}. Disable only if you know this check will be
#'   unnecessary.
#'
#' @return A list of lists of pixel \code{\link[spatstat]{im}}ages, each of which corresponds to
#'  the requested times in \code{tt}, and are named as such.\cr
#'  If \code{stob} is an object of class \code{\link{stden}}:
#'  
#'  \item{z}{
#'  Pixel images of the joint spatiotemporal density corresponding to \code{tt}.
#'  }
#'  
#'  \item{z.cond}{
#'  Pixel images of the conditional spatiotemporal density given each time in \code{tt}.
#'  }
#'  
#'  If \code{stob} is an object of class \code{\link{rrst}}:
#'  
#'  \item{rr}{
#'  Pixel images of the joint spatiotemporal relative risk corresponding to \code{tt}.
#'  }
#' 
#'  \item{rr.cond}{
#'  Pixel images of the conditional spatiotemporal relative risk given each time in \code{tt}.
#'  }
#' 
#'  \item{P}{
#'  Only present if \code{tolerate = TRUE} in the preceding call to \code{\link{spattemp.risk}}.
#'  Pixel images of the \eqn{p}-value surfaces for the joint spatiotemporal relative risk.
#'  }
#'  
#'  \item{P.cond}{
#'  Only present if \code{tolerate = TRUE} in the preceding call to \code{\link{spattemp.risk}}.
#'  Pixel images of the \eqn{p}-value surfaces for the conditional spatiotemporal relative risk.
#'  }
#'  
#' @author T.M. Davies
#'
#' @seealso \code{\link{spattemp.density}}, \code{\link{spattemp.risk}}, \code{\link{bivariate.density}}
#'
#' @references
#' Fernando, W.T.P.S. and Hazelton, M.L. (2014), Generalizing the spatial relative risk function, \emph{Spatial and Spatio-temporal Epidemiology}, \bold{8}, 1-10.
#'
#' @examples
#' 
#' \donttest{
#' data(fmd)
#' fmdcas <- fmd$cases
#' fmdcon <- fmd$controls
#' 
#' f <- spattemp.density(fmdcas,h=6,lambda=8)
#' g <- bivariate.density(fmdcon,h0=6)
#' rho <- spattemp.risk(f,g,tolerate=TRUE) 
#' 
#' f$tlim # requested slices must be in this range
#' 
#' # slicing 'stden' object
#' f.slice1 <- spattemp.slice(f,tt=50) # evaluation timestamp
#' f.slice2 <- spattemp.slice(f,tt=150.5) # interpolated timestamp
#' par(mfrow=c(2,2))
#' plot(f.slice1$z$'50')
#' plot(f.slice1$z.cond$'50')
#' plot(f.slice2$z$'150.5')
#' plot(f.slice2$z.cond$'150.5')
#' 
#' # slicing 'rrst' object
#' rho.slices <- spattemp.slice(rho,tt=c(50,150.5))
#' par(mfrow=c(2,2))
#' plot(rho.slices$rr$'50');tol.contour(rho.slices$P$'50',levels=0.05,add=TRUE)
#' plot(rho.slices$rr$'150.5');tol.contour(rho.slices$P$'150.5',levels=0.05,add=TRUE)
#' plot(rho.slices$rr.cond$'50');tol.contour(rho.slices$P.cond$'50',levels=0.05,add=TRUE)
#' plot(rho.slices$rr.cond$'150.5');tol.contour(rho.slices$P.cond$'150.5',levels=0.05,add=TRUE)
#' }
#' 
#' @export
spattemp.slice <- function(stob,tt,checkargs=TRUE){
  if(checkargs){
    if((!inherits(stob,"stden"))&&(!inherits(stob,"rrst"))) stop("'stob' must be of class \"stden\" or \"rrst\"")
    if(!is.vector(tt)||!is.numeric(tt)) stop("'tt' must be a numeric vector")
    
    if(!all(sapply(tt,function(x) x>=stob$tlim[1]) & sapply(tt,function(x) x<=stob$tlim[2]))) stop(paste("at least one requested time is outside available range of",prange(stob$tlim)))
    # was:
    # if(!inside.range(h0,aran)) stop(paste("requested 'h0' outside available range of",prange(aran)))
  }
  
  tlen <- length(tt)
  
  if(inherits(stob,"stden")){
    avail <- names(stob$z)
    z <- stob$z
    zc <- stob$z.cond
    p <- pc <- NULL
    result <- list(z=list(),z.cond=list(),P=list(),P.cond=list())
  } else {
    avail <- names(stob$rr)
    z <- stob$rr
    zc <- stob$rr.cond
    p <- stob$P
    pc <- stob$P.cond
    result <- list(rr=list(),rr.cond=list(),P=list(),P.cond=list())
  }
  
  for(i in 1:tlen){
    slc <- st.slice.single(tt[i],avail,z,zc,p,pc,warn=checkargs)
    result[[1]][[i]] <- slc$z
    result[[2]][[i]] <- slc$zc
    result[[3]][[i]] <- slc$p
    result[[4]][[i]] <- slc$pc
  }
  
  zeros <- sapply(result,length)==0
  result[which(zeros)] <- NULL
  
  names(result[[1]]) <- names(result[[2]]) <- tt
  if(length(result)>2) names(result[[3]]) <- names(result[[4]]) <- tt
  
  return(result)
}

st.slice.single <- function(V,avail,z,zc,p,pc,warn=FALSE){
  la <- length(avail)
  if(any(avail==as.character(V))){
    index <- which(avail==as.character(V))
    zres <- z[[index]]
    zcres <- zc[[index]]
    pres <- p[[index]]
    pcres <- pc[[index]]
  } else {
    marker <- as.numeric(avail)>V
    if(sum(marker)==la){
      zres <- z[[1]]
      zcres <- zc[[1]]
      pres <- p[[1]]
      pcres <- pc[[1]]
      if(warn) warning("time point lower than available range")
    } else if(sum(marker)==0){
      zres <- z[[la]]
      zcres <- zc[[la]]
      pres <- p[[la]]
      pcres <- pc[[la]]
      if(warn) warning("time point higher than available range")
    } else {
      marker <- which(marker)[1]
      mindex <- c(marker-1,marker)
      tint <- as.numeric(avail)[mindex]
      move <- (V-tint[1])/diff(tint)
      zdiff <- z[[mindex[2]]]-z[[mindex[1]]]
      zcdiff <- zc[[mindex[2]]]-zc[[mindex[1]]]
      zres <- z[[mindex[1]]]+move*zdiff
      zcres <- zc[[mindex[1]]]+move*zcdiff
      if(!is.null(p)){
        pdiff <- p[[mindex[2]]]-p[[mindex[1]]]
        pcdiff <- pc[[mindex[2]]]-pc[[mindex[1]]]
        pres <- p[[mindex[1]]]+move*pdiff
        pcres <- pc[[mindex[1]]]+move*pcdiff
      } else {
        pres <- pcres <- NULL
      }
    }
  }
  return(list(z=zres,zc=zcres,p=pres,pc=pcres))
}