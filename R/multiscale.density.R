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
  
  #' discretise window
  isrect <- is.rectangle(W)
  WM <- as.mask(W,dimyx=dimyx,xy=xy)
  insideW <- WM$m
  dimW <- WM$dim
  nr <- dimW[1]
  nc <- dimW[2]
  xstep <- WM$xstep
  ystep <- WM$ystep
  
  #' discretise spatial locations of data points
  ij <- nearest.raster.point(pp$x, pp$y, WM)
  
  #' z-mapping
  zh <- function(h,hhl) log(h)-hhl #' z map
  zhi <- function(z,hhl) exp(hhl+z) #' inverse z map
  H <- range(rep(as.vector(h.hypo.mat),2)*rep(h0fac,each=prod(dim(h.hypo.mat))),na.rm=TRUE)
  hhash.log <- log(mean(H))
  zlim <- zh(H,hhash.log)
  zrange <- c(diff(rev(zlim)),diff(zlim))
  
  #' discretise bandwidth values
  zbreaks <- seq(zrange[1], zrange[2], length=dimz+1)
  zstep <- diff(zbreaks[1:2])
  zvalues <- zbreaks[-1] - zstep/2
  kslice <- findInterval(zh(h.spec,hhash.log),zbreaks,all.inside=TRUE)
  
  #' grid coordinates (padded)
  xcol.pad <- WM$xcol[1] + xstep * (0:(2*nc-1))
  yrow.pad <- WM$yrow[1] + ystep * (0:(2*nr-1))
  z.pad <- zvalues[1] + zstep * (0:(2*dimz - 1))
  
  if(verbose) cat("Done.\nForming kernel...")
  
  #' set up kernel
  xcol.ker <- xstep * c(0:(nc-1),-(nc:1))
  yrow.ker <- ystep * c(0:(nr-1),-(nr:1))
  z.ker <- zstep * c(0:(dimz-1), -(dimz:1))
  pixarea <- xstep * ystep
  kerpixvol <- xstep * ystep * zstep
  
  #' calculating tapering values for z-coordinate
  ztap <- rep(1,2*dimz)
  if(taper){
    #' get 'zero points' to lie somewhere in the middle of the extremes beyond the raw
    #'   zrange (zrange) and within the extended zrange (zkrange)
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
  
  #' Fourier transform of kernel
  fK <- fft(Kern)
  
  if(verbose) cat("Done.\nDiscretising point locations...")
  
  rowfac <- factor(ij$row, levels=1:(2*nr))
  colfac <- factor(ij$col, levels=1:(2*nc))
  kfac <- factor(kslice, levels=1:(2*dimz))
  
  Xpad <- tapplysum(weights, list(rowfac, colfac, kfac))
  #' was:
  # Xpad <- tsum(weights, list(rowfac, colfac, kfac))
  # Xpad <- unname(unclass(Xpad))
  
  #' convolve point masses with kernel
  if(verbose) cat("Done.\nFFT of point locations...")
  fX <- fft(Xpad)
  if(verbose) cat("Inverse FFT of smoothed point locations...")
  sm <- fft(fX * fK, inverse=TRUE)/prod(dim(Xpad))
  if(verbose){
    cat("Done.\n")
    cat(paste("  [ Point convolution: maximum imaginary part=",
              signif(max(abs(Im(sm))),3), "]\n"))
  }
  
  #' identify scaling values based on range of available z-coordinates;
  #'  restrict attention to those requested by 'h0fac'
  ei <- exp(-zvalues)
  eiw <- which(ei>=h0fac[1] & ei<=h0fac[2])
  avals <- ei[eiw]
  bwvals <- avals*h0
  hlist <- list()
  for(i in 1:length(avals)) hlist[[i]] <- avals[i]*h.hypo[WM,drop=FALSE]
  
  #' turn each relevant layer of the returned FFT array into a pixel image,
  #'  giving the raw (un-edge-corrected) result for that scaling
  raw <- Re(sm)[1:nr,1:nc,1:dimz]
  rawlist <- list()
  for(i in 1:length(eiw)) rawlist[[i]] <- im(raw[,,eiw[i]],xcol=WM$xcol,yrow=WM$yrow)[W,drop=FALSE]
  names(rawlist) <- bwvals
  names(hlist) <- bwvals
  if(!intensity) rawlist <- lapply(rawlist,function(x) x/integral(x))
  result <- list(z=rev(rawlist),q=NULL,him=rev(hlist))
  WK <- elist <- NULL
  edgeW <- 1
  
  #' edge-correction
  zeta.index <- which.min((zhi(zvalues,hhash.log)-exp(hhash.log))^2) #' find slice of z corresponding to 'zero plane' (as close as possible) for placement of R^2 (works slightly better on h scale)
  if(edge=="uniform"){
    #' convolve window with kernel
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
    
    #' Uniform-type edge correction: edge weights associated with pixels
    #'  The loop below cycles through each bandwidth adjustment factor 
    #'  computed above and extracts the edge-correction surface for each.
    #'  Note that $\zeta$ is still required here since the R^2 plane won't be
    #'  *exactly* at z=0 in general due to discretisation, but will be close,
    #'  and is identified as zvalues[zeta.index]
    
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
    
    #' Construct a list of images, with names reflecting the global bandwidth scaling
    if(!intensity) lambda <- lapply(lambda,function(x) x/integral(x))
    names(elist) <- names(lambda) <- bwvals
    result$q <- rev(elist)
    result$z <- rev(lambda)
  }
  if(verbose) cat("\n")
  
  # if(taper) result$ztap <- cbind(z.ker,ztap) #' For testing/diagnostics
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