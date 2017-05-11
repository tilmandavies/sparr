edgeh <- function(bwim,pres,tres,step,W,verbose=FALSE){  
  if(pres>tres) stop("edge-correction resolution (in 'davies.baddeley[3]') exceeds spatial intensity 'resolution'")
  
  hfp <- bwim
  if(pres<tres) hfp <- as.im(interp.im,W=W,Z=bwim,dimyx=rep(pres,2))	# coarsen resolution if needed
  
  M <- as.mask(W,dimyx=rep(pres,2))
  inside <- M$m
  
  hfp[!inside] <- NA
  hfp$v[is.na(hfp$v)] <- max(hfp,na.rm=TRUE)
  hypoQ <- unique(quantile(hfp[inside],seq(0,1,step),na.rm=TRUE))
  hypocen <- hypoQ[-length(hypoQ)]+diff(hypoQ/2)
  corrQ <- as.numeric(cut(as.vector(as.matrix(hfp)),breaks=hypoQ,include.lowest=TRUE))
  if(pres<tres) corrQ[is.na(corrQ)] <- 1

  use_fftw <- fftw_available()

  lut <- matrix(FALSE, length(corrQ), length(hypoQ))
  lut[cbind(seq_along(corrQ),corrQ)] <- TRUE

  ifft_scale <- M$xstep*M$ystep/(4*pres^2)
  Mpad <- matrix(0, ncol=2*pres, nrow=2*pres)
  Mpad[1:pres, 1:pres] <- inside
  fM <- fft2d(Mpad, fftw=use_fftw)
  
  qhz <- rep(NA,pres^2)
  if(verbose) pb <- txtProgressBar(0,length(hypoQ),style=3)
  for(i in 1:length(hypocen)){
    fK <- kernel2d_fft(hypoQ[i], M$xstep, M$ystep, pres)

    con <- fft2d(fM*fK,inverse=TRUE,fftw=use_fftw)[1:pres,1:pres]
    qhz[lut[,i]] <- Mod(con[lut[,i]])*ifft_scale
    if(verbose) setTxtProgressBar(pb,i)
  }
  if(verbose) close(pb)
  
  qhz[qhz>1] <- 1
  qhz <- im(matrix(qhz,pres,pres),xcol=M$xcol,yrow=M$yrow)
  
  if(pres==tres){
    qhz[!inside] <- NA
  } else {
    qhz <- as.im(interp.im,W=W,Z=qhz,dimyx=rep(tres,2))		
  }
  
  return(qhz)	
}