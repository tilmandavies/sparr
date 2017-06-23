tol.asy.fix <- function(f,g,pooled,verbose=FALSE){

  if(verbose) cat("Initialising window...")
  if(!compatible(f$z,g$z,pooled$z)) stop("incompatible images 'f', 'g', 'pooled'... kernel estimates must be evaluated on identical domains")
  
  M <- Window(pooled$z)
  inside <- M$m
  pres <- nrow(inside)
  res2 <- 2*pres
  resseq <- 1:pres
  xcol.ker <- M$xstep*c(0:(pres-1),-rev(resseq))
  yrow.ker <- M$ystep*c(0:(pres-1),-rev(resseq))
  kerpixarea <- M$xstep*M$ystep
  len.pad <- res2^2
  Mpad <- matrix(0, ncol=2*pres, nrow=2*pres)
  Mpad[1:pres, 1:pres] <- inside
  fM <- fft(Mpad)
	
	
  fg.h0 <- mean(c(f$h0,g$h0))
  qb <- function(oo,hfac){
    hfp <- hfac*fg.h0
    densX.ker <- dnorm(xcol.ker,sd=hfp)
    densY.ker <- dnorm(yrow.ker,sd=hfp)
    Kern <- outer(densY.ker,densX.ker,"*")*kerpixarea
    con <- fft(fM*fft(Kern), inverse=TRUE)/len.pad
    qhz <- im(Mod(con[1:pres,1:pres]),xcol=oo$z$xcol,yrow=oo$z$yrow)
    qhz[qhz>1] <- 1
    qhz[!inside] <- NA
    return(as.matrix(qhz))
  }
  if(verbose) cat("Done.\nPerforming kernel*window convolution(s)...")
  
  qq <- pooled$q
  if(is.null(qq)||is.vector(qq)){
    qq <- qb(pooled,1)
  }
  qq <- as.matrix(qq)
  qq2 <- (1/(4*pi))*qb(pooled,sqrt(0.5))
  if(verbose) cat("Done.\n")
  
  RrzK <- qq2/qq^2
  denominator <- sqrt(RrzK*(1/length(f$h)+1/length(g$h)))/(fg.h0*sqrt(as.matrix(pooled$z)))
  suppressWarnings(numerator <- as.matrix(log(f$z)-log(g$z)))
  
  zstandard <- numerator/denominator
  P <- pnorm(zstandard,lower.tail=FALSE)
  return(list(num=numerator,den=denominator,zst=zstandard,p=im(P,xcol=pooled$z$xcol,yrow=pooled$z$yrow),rz=RrzK))
}