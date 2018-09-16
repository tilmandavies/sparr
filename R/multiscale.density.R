#' Multi-scale adaptive kernel density/intensity estimation
#' 
#' Computes adaptive kernel estimates of spatial density/intensity using a 3D
#' FFT for multiple global bandwidth scales.
#' 
#' Davies & Baddeley (2018) investigated computational aspects of Abramson's
#' (1982) adaptive kernel smoother for spatial (2D) data. This function is the
#' implementation of the 3D convolution via a fast-Fourier transform (FFT)
#' which allows simultaneous calculation of an adaptive kernel estimate at
#' multiple global bandwidth scales.
#' 
#' These `multiple global bandwidth scales' are computed with respect to
#' rescaling a reference value of the global bandwidth passed to the \code{h0}
#' argument. This rescaling is defined by the range provided to the argument
#' \code{h0fac}. For example, by default, the function will compute the
#' adaptive kernel estimate for a range of global bandwidths between
#' 0.25*\code{h0} and 1.5*\code{h0}. The exact numeric limits are subject to
#' discretisation, and so the returned valid range of global bandwidths will
#' differ slightly. The exact resulting range following function execution is
#' returned as the \code{h0range} element of the result, see `Value' below.
#' 
#' The distinct values of global bandwidth used (which define the
#' aforementioned \code{h0range}) and hence the total number of pixel
#' \code{\link[spatstat]{im}ages} returned depend on both the width of the span
#' \code{h0fac} and the discretisation applied to the bandwidth axis through
#' \code{dimz}. Increasing this z-resolution will provide more pixel images and
#' hence greater numeric precision, but increases computational cost. The
#' returned pixel \code{\link[spatstat]{im}ages} that represent the multiscale
#' estimates are stored in a named list (see `Value'), whose names reflect the
#' corresponding distinct global bandwidth. See `Examples' for the easy way to
#' extract these distinct global bandwidths.
#' 
#' The user can request an interpolated density/intensity estimate for any
#' global bandwidth value within \code{h0range} by using the
#' \code{\link{multiscale.slice}} function, which returns an object of class
#' \code{\link{bivden}}.
#' 
#' @aliases multiscale.density msden
#' 
#' @param pp An object of class \code{\link[spatstat]{ppp}} giving the observed
#'   2D data set to be smoothed.
#' @param h0 Reference global bandwidth for adaptive smoothing; numeric value >
#'   0. Multiscale estimates will be computed by rescaling this value as per
#'   \code{h0fac}.
#' @param hp Pilot bandwidth (scalar, numeric > 0) to be used for fixed
#'   bandwidth estimation of the pilot density. If \code{NULL} (default), it will
#'   take on the value of \code{h0}. Ignored when \code{pilot.density} is
#'   supplied as a pre-defined pixel image.
#' @param h0fac A numeric vector of length 2 stipulating the span of the global
#'   bandwidths in the multiscale estimates. Interpreted as a multiplicative
#'   factor on \code{h0}. See `Details'.
#' @param edge Character string dictating edge correction. \code{"uniform"}
#'   (default) corrects based on evaluation grid coordinate. Setting \code{edge="none"}
#'   requests no edge correction.
#' @param resolution Numeric value > 0. Resolution of evaluation grid in the
#'   spatial domain; the densities/intensities will be returned on a
#'   [\code{resolution} \eqn{\times}{x} \code{resolution}] grid.
#' @param dimz Resolution of z- (rescaled bandwidth)-axis in the trivariate
#'   convolution. Higher values increase precision of the multiscale estimates at
#'   a computational cost. See `Details'.
#' @param gamma.scale Scalar, numeric value > 0; controls rescaling of the
#'   variable bandwidths. Defaults to the geometric mean of the bandwidth factors
#'   given the pilot density (as per Silverman, 1986). See the documentation for
#'   \code{\link{bivariate.density}}.
#' @param trim Numeric value > 0; controls bandwidth truncation for adaptive
#'   estimation. See the documentation for \code{\link{bivariate.density}}.
#' @param intensity Logical value indicating whether to return an intensity
#'   estimate (integrates to the sample size over the study region), or a density
#'   estimate (default, integrates to 1).
#' @param pilot.density An optional pixel image (class
#'   \code{\link[spatstat]{im}}) giving the pilot density to be used for
#'   calculation of the variable bandwidths in adaptive estimation, \bold{or} a
#'   \code{\link[spatstat]{ppp.object}} giving the data upon which to base a
#'   fixed-bandwidth pilot estimate using \code{hp}. See the documentation for
#'   \code{\link{bivariate.density}}.
#' @param xy Optional alternative specification of the spatial evaluation grid;
#'   matches the argument of the same tag in \code{\link[spatstat]{as.mask}}. If
#'   supplied, \code{resolution} is ignored.
#' @param taper Logical value indicating whether to taper off the trivariate
#'   kernel outside the range of \code{h0*h0fac} in the scale space; see Davies &
#'   Baddeley (2018). Keep at the default \code{TRUE} if you don't know what this
#'   means.
#' @param verbose Logical value indicating whether to print function progress.
#' 
#' @return An object of class \code{"msden"}. This is very similar to a
#' \code{\link{bivden}} object, with lists of pixel
#' \code{\link[spatstat]{im}}ages in the \code{z}, \code{him}, and \code{q}
#' components (instead of standalone images).
#' \item{z}{A list of the resulting
#'   density/intensity estimates; each member being a pixel image object of class
#'   \code{\link[spatstat]{im}}. They are placed in increasing order of the
#'   discretised values of \code{h0}.}
#' \item{h0}{A copy of the reference value of \code{h0} used.}
#' \item{h0range}{A vector of length 2 giving the actual range
#'   of global bandwidth values available (inclusive).}
#' \item{hp}{A copy of the value of \code{hp} used.}
#' \item{h}{A numeric vector of length equal to the
#'   number of data points, giving the bandwidth used for the corresponding
#'   observation in \code{pp} with respect to the reference global bandwidth
#'   \code{h0}.}
#' \item{him}{A list of pixel images (class \code{\link[spatstat]{im}}),
#'   corresponding to \code{z}, giving the
#'   `hypothetical' Abramson bandwidth at each pixel coordinate conditional upon
#'   the observed data and the global bandwidth used.}
#' \item{q}{Edge-correction weights; list of pixel \code{\link[spatstat]{im}}ages
#'   corresponding to \code{z} if \code{edge = "uniform"}, and \code{NULL} if
#'   \code{edge = "none"}.}
#' \item{gamma}{The numeric value of \code{gamma.scale} used in scaling the bandwidths.}
#' \item{geometric}{The geometric mean of the
#'   untrimmed variable bandwidth factors. This will be identical to \code{gamma}
#'   if \code{gamma.scale = "geometric"} as per default.}
#' \item{pp}{A copy of the \code{\link[spatstat]{ppp.object}} initially passed to the
#'   \code{pp} argument, containing the data that were smoothed.}
#' 
#' @author T.M. Davies and A. Baddeley
#' 
#' @seealso \code{\link{bivariate.density}}, \code{\link{multiscale.slice}}
#'
#' @references
#' 
#' Abramson, I. (1982). On bandwidth variation in kernel estimates
#' --- a square root law, \emph{Annals of Statistics}, \bold{10}(4),
#' 1217-1223.
#' 
#' Davies, T.M. and Baddeley A. (2018), Fast computation of
#' spatially adaptive kernel estimates, \emph{Statistics and Computing}, \bold{28}(4), 937-956.
#' 
#' Silverman, B.W. (1986), \emph{Density Estimation for Statistics and Data Analysis},
#' Chapman & Hall, New York.
#' 
#' @examples
#' \donttest{
#' data(chorley) # Chorley-Ribble data (package 'spatstat')
#' ch.multi <- multiscale.density(chorley,h0=1)
#' plot(ch.multi)
#' 
#' ch.pilot <- bivariate.density(chorley,h0=0.75) # with pre-defined pilot density
#' ch.multi2 <- multiscale.density(chorley,h0=1,pilot.density=ch.pilot$z)
#' plot(ch.multi2)
#' 
#' data(pbc)
#' # widen h0 scale, increase z-axis resolution
#' pbc.multi <- multiscale.density(pbc,h0=2,hp=1,h0fac=c(0.25,2.5),dimz=128) 
#' plot(pbc.multi)
#' }
#' @export
multiscale.density <- function(pp,h0,hp=NULL,h0fac=c(0.25,1.5),edge=c("uniform","none"),resolution=128,dimz=64,gamma.scale="geometric",trim=5,intensity=FALSE,pilot.density=NULL,xy=NULL,taper=TRUE,verbose=TRUE){
  if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")
  
  if(verbose) cat("Initialising...")
  
  n <- npoints(pp)
  W <- Window(pp)
  if(!is.null(xy)){
    xy <- checkxy(xy)
    dimyx <- NULL
    resolution <- length(xy$x)
  } else {
    resolution <- checkit(resolution,"'resolution'")
    dimyx <- rep(resolution,2)
  }
  dimz <- checkit(dimz,"'dimz'")
  h0fac <- checkran(h0fac,"h0fac")
  h0 <- checkit(h0,"'h0'")

  if(is.null(hp)) hp <- h0
  else hp <- checkit(hp,"'hp'")
  
  edge <- checkedge(edge,v=0)
  
  if(!is.null(trim)&&!is.na(trim)) trim <- checkit(trim,"'trim'")
  
  pd <- pilot.density
  pilot.data <- pp
  if(!is.null(pd)){
    if(is.im(pd)){
      if(is.null(xy)){
        if(!all(dim(pd)==resolution)) stop("'pilot.density' image resolution must strictly have 'resolution' x 'resolution' pixels")
      } else {
        if((!all(pd$xcol==xy$x))||(!all(pd$yrow==xy$y))) stop("'pilot.density' xcol and yrow must strictly match coords in 'xy'")
      }
      pilot.density[pd<=0] <- min(pd[pd>0])
    } else if(is.ppp(pd)){
      pilot.data <- pd
      if(!identical_windows(W,Window(pilot.data))) stop("'pilot.density' window must be identical to 'pp' window")
      pilot.density <- density(pilot.data,sigma=hp,edge=(edge=="uniform"),diggle=FALSE,dimyx=dimyx,xy=xy,positive=TRUE)
    } else {
      stop("'pilot.density' must be an object of class \"im\" or \"ppp\"")
    }
  } else {
    pilot.density <- density(pp,sigma=hp,edge=(edge=="uniform"),diggle=FALSE,dimyx=dimyx,xy=xy,positive=TRUE)
  }
  
  pilot.density.spec <- safelookup(pilot.density,pp,warn=FALSE)
  pi.int <- integral(pilot.density)
  pilot.density <- pilot.density/pi.int
  pilot.density.spec <- pilot.density.spec/pi.int
  pspec <- pilot.density.spec^(-0.5)
  gamma <- processgamma(gamma.scale,safelookup(pilot.density,pilot.data,warn=FALSE))
  gs <- gspd <- exp(mean(log(pspec))) 
  if(!is.null(pd)) gspd <- exp(mean(log(safelookup(pilot.density,pilot.data,warn=FALSE)^(-0.5))))
  
  h.spec <- h0*pmin(pspec,trim*gspd)/gamma
  h.hypo <- h0*im(matrix(pmin(as.vector(as.matrix(pilot.density^(-0.5))),trim*gspd),resolution,resolution)/gamma,xcol=pilot.density$xcol,yrow=pilot.density$yrow)
  h.hypo.mat <- as.matrix(h.hypo)
  
  marks(pp) <- h.spec
  weights <- rep(1,n)
  
  if(verbose) cat("Done.\nDiscretising...")
  
  # discretise window
  isrect <- is.rectangle(W)
  WM <- as.mask(W,dimyx=dimyx,xy=xy)
  insideW <- WM$m
  dimW <- WM$dim
  nr <- dimW[1]
  nc <- dimW[2]
  xstep <- WM$xstep
  ystep <- WM$ystep
  
  # discretise spatial locations of data points
  ij <- nearest.raster.point(pp$x, pp$y, WM)
  
  # z-mapping
  zh <- function(h,hhl) log(h)-hhl #' z map
  zhi <- function(z,hhl) exp(hhl+z) #' inverse z map
  H <- range(rep(as.vector(h.hypo.mat),2)*rep(h0fac,each=prod(dim(h.hypo.mat))),na.rm=TRUE)
  hhash.log <- log(mean(H))
  zlim <- zh(H,hhash.log)
  zrange <- c(diff(rev(zlim)),diff(zlim))
  
  # discretise bandwidth values
  zbreaks <- seq(zrange[1], zrange[2], length=dimz+1)
  zstep <- diff(zbreaks[1:2])
  zvalues <- zbreaks[-1] - zstep/2
  kslice <- findInterval(zh(h.spec,hhash.log),zbreaks,all.inside=TRUE)
  
  # grid coordinates (padded)
  xcol.pad <- WM$xcol[1] + xstep * (0:(2*nc-1))
  yrow.pad <- WM$yrow[1] + ystep * (0:(2*nr-1))
  z.pad <- zvalues[1] + zstep * (0:(2*dimz - 1))
  
  if(verbose) cat("Done.\nForming kernel...")
  
  # set up kernel
  xcol.ker <- xstep * c(0:(nc-1),-(nc:1))
  yrow.ker <- ystep * c(0:(nr-1),-(nr:1))
  z.ker <- zstep * c(0:(dimz-1), -(dimz:1))
  pixarea <- xstep * ystep
  kerpixvol <- xstep * ystep * zstep
  
  # calculating tapering values for z-coordinate
  ztap <- rep(1,2*dimz)
  if(taper){
    # get 'zero points' to lie somewhere in the middle of the extremes beyond the raw
    #   zrange (zrange) and within the extended zrange (zkrange)
    zkrange <- range(z.ker)
    zlo <- zkrange[1]+(zrange[1]-zkrange[1])/2
    zup <- zkrange[2]-(zkrange[2]-zrange[2])/2
    ztap <- pmin(taperoff(z.ker,zlo,zrange[1],"cosine"),taperoff(z.ker,zup,zrange[2],"cosine"))
  }
  
  weights <- weights/pixarea
  Kern <- array(0, dim=2*c(nr, nc, dimz))
  hk <- zhi(-z.ker,hhash.log)
  for(k in 1:(2*dimz)) {
    densX.ker <- dnorm(xcol.ker, sd=hk[k])
    densY.ker <- dnorm(yrow.ker, sd=hk[k])
    Kern[,,k] <- outer(densY.ker, densX.ker, "*") * pixarea * ztap[k] #' z-tapering inserted here
  }
  
  if(verbose) cat("Done.\nTaking FFT of kernel...")
  
  # Fourier transform of kernel
  fK <- fft(Kern)
  
  if(verbose) cat("Done.\nDiscretising point locations...")
  
  rowfac <- factor(ij$row, levels=1:(2*nr))
  colfac <- factor(ij$col, levels=1:(2*nc))
  kfac <- factor(kslice, levels=1:(2*dimz))
  
  Xpad <- tapplysum(weights, list(rowfac, colfac, kfac))
  # was:
  # Xpad <- tsum(weights, list(rowfac, colfac, kfac))
  # Xpad <- unname(unclass(Xpad))
  
  # convolve point masses with kernel
  if(verbose) cat("Done.\nFFT of point locations...")
  fX <- fft(Xpad)
  if(verbose) cat("Inverse FFT of smoothed point locations...")
  sm <- fft(fX * fK, inverse=TRUE)/prod(dim(Xpad))
  if(verbose){
    cat("Done.\n")
    cat(paste("  [ Point convolution: maximum imaginary part=",
              signif(max(abs(Im(sm))),3), "]\n"))
  }
  
  # identify scaling values based on range of available z-coordinates;
  #  restrict attention to those requested by 'h0fac'
  ei <- exp(-zvalues)
  eiw <- which(ei>=h0fac[1] & ei<=h0fac[2])
  if(length(eiw)==0) stop("'h0fac' too narrow for 'dimz' -- increase either one or both")
  avals <- ei[eiw]
  bwvals <- avals*h0
  hlist <- list()
  for(i in 1:length(avals)) hlist[[i]] <- avals[i]*h.hypo[WM,drop=FALSE]
  
  
  
  
  # turn each relevant layer of the returned FFT array into a pixel image,
  #  giving the raw (un-edge-corrected) result for that scaling
  raw <- Re(sm)[1:nr,1:nc,1:dimz]
  rawlist <- list()
  for(i in 1:length(eiw)) rawlist[[i]] <- im(raw[,,eiw[i]],xcol=WM$xcol,yrow=WM$yrow)[W,drop=FALSE]
  names(rawlist) <- bwvals
  names(hlist) <- bwvals
  if(!intensity) rawlist <- lapply(rawlist,function(x) x/integral(x))
  result <- list(z=rev(rawlist),q=NULL,him=rev(hlist))
  WK <- elist <- NULL
  edgeW <- 1
  
  # edge-correction
  zeta.index <- which.min((zhi(zvalues,hhash.log)-exp(hhash.log))^2) #' find slice of z corresponding to 'zero plane' (as close as possible) for placement of R^2 (works slightly better on h scale)
  if(edge=="uniform"){
    # convolve window with kernel
    Wpad <- array(0, dim=2*c(nr, nc, dimz))
    Wpad[1:nr, 1:nc, zeta.index] <- WM$m * pixarea #' insert W plane here
    if(verbose) cat("FFT of window...")
    fW <- fft(Wpad)
    if(verbose) cat("Inverse FFT of smoothed window...")
    WK <- fft(fW * fK, inverse=TRUE)/prod(dim(Wpad))
    if(verbose){
      cat("Done.\n")
      cat(paste("  [ Window convolution: maximum imaginary part=",
                signif(max(abs(Im(WK))),3), "]\n"))
      cat("Looking up edge correction weights...\n")
    }
    
    elist <- lambda <- list()   
    
    # Uniform-type edge correction: edge weights associated with pixels
    #  The loop below cycles through each bandwidth adjustment factor 
    #  computed above and extracts the edge-correction surface for each.
    #  Note that $\zeta$ is still required here since the R^2 plane won't be
    #  *exactly* at z=0 in general due to discretisation, but will be close,
    #  and is identified as zvalues[zeta.index]
    
    for(i in 1:length(avals)){
      if(verbose) cat(paste(i," "))
      
      # adj <- avals[i]
      edgeW <- bwW <- hlist[[i]]#adj*h.hypo[WM, drop=FALSE]
      
      kim <- eval.im(findInterval(zvalues[zeta.index]+hhash.log-log(bwW),zbreaks,all.inside=TRUE))
      iim <- as.im(row(bwW),W=bwW)
      jim <- as.im(col(bwW),W=bwW)
      df <- pairs(iim,jim,kim,plot=FALSE)
      
      edgeW[] <- Re(WK)[as.matrix(df)]/pixarea
      elist[[i]] <- edgeW
      lambda[[i]] <- rawlist[[i]]/edgeW
    }
    
    # Construct a list of images, with names reflecting the global bandwidth scaling
    if(!intensity) lambda <- lapply(lambda,function(x) x/integral(x))
    names(elist) <- names(lambda) <- bwvals
    result$q <- rev(elist)
    result$z <- rev(lambda)
  }
  if(verbose) cat("\n")
  
  # if(taper) result$ztap <- cbind(z.ker,ztap) # For testing/diagnostics
  result$h <- h.spec
  result$h0range <- range(bwvals)
  result$gamma <- gamma
  result$geometric <- gs
  result$pp <- pp
  result$hp <- hp
  result$h0 <- h0
  
  result <- result[c("z","h0","h0range","hp","h","him","q","gamma","geometric","pp")]
  class(result) <- "msden"
  return(result)
}
