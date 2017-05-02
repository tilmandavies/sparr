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
  imlist <- lapply(hu,subpix,plocs=nearest.raster.point(x$x,x$y,w=WM),hc=hc,WM=WM,weights=weights)
  xcol.pad <- WM$xcol[1]+WM$xstep*(0:(res2-1))
  yrow.pad <- WM$yrow[1]+WM$ystep*(0:(res2-1))
  xcol.ker <- WM$xstep*c(0:(resolution-1),-rev(resseq))
  yrow.ker <- WM$ystep*c(0:(resolution-1),-rev(resseq))
  kerpixarea <- WM$xstep*WM$ystep
  len.pad <- res2^2
  resultlist <- list()
  result <- im(matrix(0,resolution,resolution),xcol=WM$xcol,yrow=WM$yrow)
  if(verbose) pb <- txtProgressBar(0,U)
  for(i in 1:U){
    nn <- sum(hc==hu[i])
    Xi <- imlist[[i]]
    xi <- as.matrix(Xi)
    
    zpad <- matrix(0,res2,res2)
    zpad[resseq,resseq] <- xi
    densX.ker <- dnorm(xcol.ker,sd=hu[i])
    densY.ker <- dnorm(yrow.ker,sd=hu[i])
    Kern <- outer(densY.ker,densX.ker,"*")*kerpixarea
    
    fK <- fft(Kern)
    fZ <- fft(zpad)
    sm <- fft(fZ*fK,inverse=TRUE)/len.pad
    smo <- im(Re(sm)[resseq,resseq],xcol.pad[resseq],yrow.pad[resseq])
    smo$v <- smo$v/kerpixarea
    
    resultlist[[i]] <- smo
    result <- result + smo 
    if(verbose) setTxtProgressBar(pb,i)
  }
  if(verbose) close(pb)
  result[!insideW] <- NA
  
  if(edge && !diggle) result <- result/edg
  if(!intensity) result <- result/integral(result)
  
  return(list(result=result,rlist=resultlist,edg=edg))
}