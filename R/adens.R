point_image_by_bw <- function(bw_cat, bw, points, weights, WM) {
  # bw is the binned bandwidths evaluated at the points
  # bw_cat is the unique, sorted list of bins

  # We want to create a surface per bin, with each cell containing the
  # (weighted) sum of points that use that bin.

  # We can do this efficiently by iterating over the points once, and just
  # adding each point to the appropriate surface.

  # Create a map from the bandwidths for each point to our bandwidth categories
  bw_map <- match(bw, bw_cat)

  # Create output matrices for each surface. We create them twice as large
  # so they are padded all ready for the FFT that is coming later
  im_bw <- vector('list', length(bw_cat))
  for (i in 1:length(bw_cat)) {
    im_bw[[i]] <- matrix(0, WM$dim[1]*2, WM$dim[2]*2)
  }

  # Iterate over the points, using bw_map to map point to surface, and fill them in
  for (i in 1:length(bw_map)) {
    cat = bw_map[i]
    x = points$row[i]
    y = points$col[i]
    im_bw[[cat]][x,y] = im_bw[[cat]][x,y] + weights[i]
  }

  im_bw
}

adens <- function(x,bwim,bwpts,resolution,edge,diggle,weights,intensity,hstep,qstep,qres,verbose){
  n <- npoints(x)
  hc <- gethcats(bwpts,step=hstep)
  hu <- sort(unique(hc))
  U <- length(hu)
  
  W <- Window(x)
  WM <- as.mask(W,dimyx=rep(resolution,2))
  insideW <- WM$m
  res2 <- 2*resolution
  resseq <- 1:resolution
  
  edg <- im(matrix(1,resolution,resolution),xcol=WM$xcol,yrow=WM$yrow)
  if(edge){
  	edg <- edgeh(bwim,pres=qres,tres=resolution,step=qstep,W=W)
  }
  
  if(is.null(weights)) weights <- rep(1,n)
  
  if(edge&&diggle) digw <- 1/safelookup(edg,x,warn=FALSE)
  else digw <- rep(1,n)
  
  weights <- weights*digw
  imlist <- point_image_by_bw(hu, hc, nearest.raster.point(x$x,x$y,w=WM), weights, WM)

  xcol.pad <- WM$xcol[1]+WM$xstep*(0:(res2-1))
  yrow.pad <- WM$yrow[1]+WM$ystep*(0:(res2-1))
  xcol.ker <- WM$xstep*c(0:(resolution-1),-rev(resseq))
  yrow.ker <- WM$ystep*c(0:(resolution-1),-rev(resseq))
  len.pad <- res2^2
  resultlist <- list()
  result <- matrix(0,resolution,resolution)
  if(verbose) pb <- txtProgressBar(0,U)
  for(i in 1:U){
    densX.ker <- dnorm(xcol.ker,sd=hu[i])
    densY.ker <- dnorm(yrow.ker,sd=hu[i])

    fK <- Re(fft(densY.ker)) %*% t(Re(fft(densX.ker)))
    fZ <- fft(imlist[[i]])
    sm <- fft(fZ*fK,inverse=TRUE)/len.pad
    smo <- Re(sm[resseq,resseq])

    resultlist[[i]] <- im(smo,xcol.pad[resseq],yrow.pad[resseq])
    result <- result + smo 
    if(verbose) setTxtProgressBar(pb,i)
  }
  if(verbose) close(pb)
  result[!insideW] <- NA
  result <- im(result,xcol=WM$xcol,yrow=WM$yrow)

  if(edge && !diggle) result <- result/edg
  if(!intensity) result <- result/integral(result)
  
  return(list(result=result,rlist=resultlist,edg=edg))
}