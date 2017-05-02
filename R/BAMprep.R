
BAMprep <- function(cases,controls,lambda,erode,res){
  W <- Window(cases)
  WM <- as.mask(W,dimyx=res)
  yx <- expand.grid(WM$yrow,WM$xcol)
  
  if(is.na(erode)){
    xyin <- inside.owin(x=yx[,2],y=yx[,1],w=W)
  } else {
    ewin <- erosion(W,erode*lambda,polygonal=TRUE)
    xyin <- inside.owin(x=yx[,2],y=yx[,1],w=ewin)
  }
  
  fdd <- nu.dashdash(cbind(cases$x,cases$y),yx[,2],yx[,1],lambda,xyin)
  gdd <- nu.dashdash(cbind(controls$x,controls$y),yx[,2],yx[,1],lambda,xyin)
  
  inside <- WM$m
  res2 <- 2*res
  resseq <- 1:res
  xcol.ker <- WM$xstep*c(0:(res-1),-rev(resseq))
  yrow.ker <- WM$ystep*c(0:(res-1),-rev(resseq))
  kerpixarea <- WM$xstep*WM$ystep
  len.pad <- res2^2
  Mpad <- matrix(0, ncol=2*res, nrow=2*res)
  Mpad[1:res,1:res] <- inside
  fM <- fft(Mpad)
  
  return(list(cas=cases,con=controls,xyin=xyin,fdd=fdd,gdd=gdd,res=res,ero=erode,fM=fM,xk=xcol.ker,yk=yrow.ker,kpa=kerpixarea,lp=len.pad,WM=WM))
}

nu.dashdash <- function(data,xs,ys,lambda,xyin){
  res2 <- length(xs)
  n <- nrow(data)
  result <- rep(NA,res2)
  for(i in 1:res2){
    if(!xyin[i]) next
    xyir <- cbind(rep(xs[i],n),rep(ys[i],n))
    data.xydiff2 <- rowSums((data-xyir)^2)
    kdashdash <- (data.xydiff2/(2*lambda)-1)*exp(data.xydiff2/(-2*lambda^2))
    result[i] <- (1/(n*pi*lambda^4))*sum(kdashdash)
  }
  return(result)
}