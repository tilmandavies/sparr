#' @export
BOOT.spattemp <- function(pp,tt=NULL,tlim=NULL,eta=NULL,nu=NULL,
                          sedge=c("uniform","none"),tedge=sedge,
                          ref.density=NULL,sres=64,tres=sres,
                          start=NULL,verbose=TRUE){

  if(!inherits(pp,"ppp")) stop("'pp' must be of spatstat class \"ppp\"; see ?ppp")
  
  if(verbose) cat("Initialising...")
  
  n <- npoints(pp)
  W <- Window(pp)
  WM <- as.mask(W,dimyx=rep(sres,2))
  
  if(is.null(tt)) tt <- marks(pp)
  tt <- checktt(tt)
  if(length(tt)!=n) stop(paste("Length of temporal vector does not match number of spatial observations\n   npoints(pp) = ",n,"; length(tt) = ",length(tt),sep=""))
  if(is.null(tlim)) tlim <- range(tt)
  tlim <- checkranin(tlim,tt,"tlim")
  
  if(is.null(tres)){
    tcw <- 1
    kt <- tlim <- c(floor(tlim[1]),ceiling(tlim[2]))
    grt <- tlim[1]:tlim[2]
    tres <- length(grt)
  } else {
    tres <- checkit(tres,"'tres'")
    tcw <- diff(tlim)/tres
    grt <- tlim[1]+0.5*tcw+(0:(tres-1))*tcw
    kt <- c(tlim[1]+0.5*tcw,tlim[2]-0.5*tcw)
  }
  
  # LIMIT OPTIMISATION? CURRENTLY UNIMPLEMENTED.
  # # set h limits if unsupplied
  # if(is.null(hlim)){
  #   ppu <- pp
  #   marks(ppu) <- NULL
  #   md <- min(nndist(unique(ppu)))
  #   hlim <- c(md,max(md*50,min(diff(W$xrange),diff(W$yrange))/6))
  # } else {
  #   hlim <- checkran(hlim,"'hlim'")
  # }
  # 
  # # set lambda limits if unsupplied
  # if(is.null(lambdalim)){
  #   ttu <- unique(tt)
  #   ttd <- outer(ttu,ttu,"-")
  #   mt <- min(ttd[lower.tri(ttd)])
  #   lamlim <- c(mt,max(mt*50,diff(tlim)))
  # } else {
  #   lamlim <- checkran(lambdalim,"'lamlim'")
  # }
  
  sedg <- checkedge(sedge,v=0)=="uniform"
  tedg <- checkedge(tedge,v=0)=="uniform"
  
  # FOR FUTURE IMPLEMENTATION
  # if(!is.na(parallelise)){
  #   if(!is.numeric(parallelise)) stop("'parallelise' must be numeric")
  #   if(is.null(parallelise)) parallelise <- NA
  #   parallelise <- round(parallelise[1])
  # }
  
  inside <- WM$m
  inn <- which(as.vector(inside))
  evalyx <- as.matrix(expand.grid(WM$yrow,WM$xcol))
  evalyx.redu <- evalyx[inn,]
  GN <- length(inn)
  
  if(is.null(ref.density)){
    
    if(is.null(eta)) eta <- OS(pp)
    else eta <- checkit(eta,"'eta'")
    
    if(is.null(nu)) nu <- ((243/(2*sqrt(pi)))/(35*n))^(1/5)*min(sd(tt),IQR(tt)/1.34) # Univariate oversmoothing bandwidth
    else nu <- checkit(nu,"'nu'")
    
    # if(sedg){
    #   ref.density <- spattemp.density(pp,h=eta,tt=tt,lambda=nu,tlim=tlim,sedge="uniform",)
    #
    #   d.eta <- density(pp,eta,edge=TRUE,positive=TRUE,dimyx=sres,spill=1)
    #   d.eta$edg[d.eta$edg>1] <- 1
    #   d.etadens <- d.eta$raw/d.eta$edg
    #   d.etaint <- integral(d.etadens)
    #   d.etadens <- as.matrix(d.etadens)[inn]/d.etaint
    #   epsilon.eta <- safelookup(d.eta$edg,pp,warn=FALSE)
    # } else {
    #   # d.etadens <- as.matrix(d.eta$raw/integral(d.eta$raw))[inn]
    #   epsilon.eta <- rep(1,n)
    # }
    # 
    # if(tedg) epsilon.nu <- pnorm(tlim[2],mean=tt,sd=nu) - pnorm(tlim[1],mean=tt,sd=nu)
    # else epsilon.nu <- rep(1,n)

    ref.density <- spattemp.density(pp,h=eta,tt=tt,lambda=nu,tlim=tlim,sedge=sedge,tedge=tedge,sres=sres,tres=tres,verbose=FALSE)
    
        
  } else {
    if(!inherits(ref.density,"stden")) stop("'ref.density' must be of class \"stden\"; see ?spattemp.density")
    if(!compatible(as.im(WM),ref.density$z[[1]])) stop("'ref.density' must be evaluated on identical spatial domain as 'Window(pp)' given 'sres'")
    if(sedg&&is.null(ref.density$qs)) stop("'ref.density' spatial edge-correction must exist if sedge = \"uniform\"")
    
    # spatial matching
    eta <- ref.density$h
    nu <- ref.density$lambda
    
    # temporal matching
    # reftres <- length(ref.density$grt)
    # if(tres!=reftres) stop("'ref.density' temporal resolution (currently ",reftres,") must match 'tres' (currently ",tres,")",sep="")
    # if(!all(ref.density$tlim==range(grt))) stop("ref.density$tlim must be identical to 'tlim'")
    # if(tedg&&is.null(ref.density$qt)) stop("'ref.density' temporal edge-correction must exist if tedge = \"uniform\"")
  }
  
  if(!is.null(ref.density$qs)){
    epsilon.eta <- safelookup(ref.density$qs,pp,warn=FALSE)
    use_fftw <- fftw_available()
    ifft_scale <- WM$xstep*WM$ystep/(4*sres^2)
    Mpad <- matrix(0,2*sres,2*sres)
    Mpad[1:sres,1:sres] <- inside
    fM <- fft2d(Mpad,fftw=use_fftw)
  } else {
    epsilon.eta <- rep(1,n)
    ifft_scale <- use_fftw <- fM <- NA
  }
  
  if(!is.null(ref.density$qt)) epsilon.nu <- pnorm(tlim[2],mean=tt,sd=nu) - pnorm(tlim[1],mean=tt,sd=nu)
  else epsilon.nu <- rep(1,n)

  xs <- ref.density$z[[1]]$xstep
  ys <- ref.density$z[[1]]$ystep
  ts <- ref.density$tgrid[2]-ref.density$tgrid[1]
  sqz <- lapply(ref.density$z,function(x) x^2)
  boot3 <- sum(Reduce("+",sqz)*xs*ys*ts,na.rm=TRUE)
  
  if(verbose) cat("Done.\nOptimising...\n")
  
  if(is.null(start)){
    start <- c(eta,nu)
  } else {
    if(any(start<=0)) stop("Invalid starting values in 'start'")
  }
  
  result <- optim(par=start,boot.opt.spattemp.fix,sedg=sedg,tedg=tedg,WM=WM,tlim=tlim,
                  sres=sres,tres=tres,fM=fM,ifft_scale=ifft_scale,inn=inn,GN=GN,
                  evalyx.redu=evalyx.redu,evalt=grt,pp=pp,tt=tt,epsilon.eta=epsilon.eta,
                  epsilon.nu=epsilon.nu,eta=eta,nu=nu,n=n,boot3=boot3,use_fftw=use_fftw,
                  parallelise=NA,verbose=verbose)$par
  
  if(verbose) cat("Done.\n")
  
  names(result) <- c("h","lambda")
  return(result)
}


boot.opt.spattemp.fix <- function(hlam,sedg,tedg,WM,tlim,sres,tres,fM,ifft_scale,inn,GN,evalyx.redu,evalt,
                                  pp,tt,epsilon.eta,epsilon.nu,eta,nu,n,boot3,use_fftw,parallelise,verbose){
  if(any(hlam<=0)) return(NA)
  h <- hlam[1]
  lam <- hlam[2]
  
  if(verbose) cat("h =",h,"\b; lambda =",lam,"\n")
  
  if(sedg){
    fK.h <- kernel2d_fft(h,WM$xstep,WM$ystep,sres)
    fK.con <- fft2d(fM*fK.h,inverse=TRUE,use_fftw)[1:sres,1:sres]
    edg.h <- Mod(fK.con)*ifft_scale
    edg.h[edg.h>1] <- 1
    GE.h <- edg.h[inn]
  } else {
    GE.h <- rep(1,length(inn))
  }
  
  if(tedg) GE.lam <- pnorm(tlim[2],mean=evalt,sd=lam) - pnorm(tlim[1],mean=evalt,sd=lam)
  else GE.lam <- rep(1,tres)
  
  if(is.na(parallelise)){
    bs2temp <- bs1 <- bs2 <- matrix(NA,n,GN)
    for(i in 1:GN){
      evx <- rep(evalyx.redu[i,2],n)-pp$x
      evy <- rep(evalyx.redu[i,1],n)-pp$y
      bs2temp[,i] <- epsilon.eta^(-1)*kernel2d(evx,evy,sqrt(h^2+eta^2))
      bs2[,i] <- epsilon.eta^(-1)*kernel2d(evx,evy,eta)
      bs1[,i] <- epsilon.eta^(-1)*kernel2d(evx,evy,sqrt(0.5*h^2+eta^2))
    }

    bt2temp <- bt1 <- bt2 <- matrix(NA,n,tres)
    b1c <- b2c <- matrix(NA,GN,tres)
    for(j in 1:tres){
      et <- rep(evalt[j],n)-tt
      bt2temp[,j] <- epsilon.nu^(-1)*dnorm(et,0,sqrt(lam^2+nu^2))
      bt2[,j] <- epsilon.nu^(-1)*dnorm(et,0,nu)
      bt1[,j] <-  epsilon.nu^(-1)*dnorm(et,0,sqrt(0.5*lam^2+nu^2))
    
      b2temp.x.j <- colSums(bs2temp*matrix(rep(bt2temp[,j],GN),n,GN))
      b2c[,j] <- b2temp.x.j*colSums(bs2*matrix(rep(bt2[,j],GN),n,GN))
      b1c[,j] <- (1/(8*pi^1.5*h^2*lam))*colSums(bs1*matrix(rep(bt1[,j],GN),n,GN)) + ((n-1)/n)*b2temp.x.j^2
    }
  } else {
    stop("'parallelise' unimplemented")
  }
  
  cella <- WM$xstep*WM$ystep*(evalt[2]-evalt[1])

  dexi <- matrix(rep(GE.h,tres),GN,tres)
  deti <- matrix(rep(GE.lam,GN),GN,tres,byrow=TRUE)
  
  boot1 <- deti^(-2)*dexi^(-2)*b1c
  boot2 <- -2*deti^(-1)*dexi^(-1)*b2c

  return(sum(n^(-2)*(boot1+boot2)*cella,na.rm=TRUE)+boot3)
}




# boot.hlam <- function(bands,ppt,g=bands[1],gam=bands[2],tlim=NULL,sres=64,tres=sres,parallelise=NA){
#   if(any(bands<=0)) return(NA)
#   
#   W <- Window(ppt)
#   Wm <- as.mask(W,dimyx=rep(sres,2))
#   n <- npoints(ppt)
#   ppdat <- cbind(ppt$x,ppt$y)
#   tt <- marks(ppt)
#   if(is.null(tlim)) tlim <- range(tt)
#   h <- bands[1]
#   lam <- bands[2]
#   
#   dh <- density(ppt,h,diggle=FALSE,positive=TRUE,dimyx=c(sres,sres),spill=1)
#   dh$edg[dh$edg>1] <- 1
#   dg <- density(ppt,g,diggle=TRUE,positive=TRUE,dimyx=c(sres,sres),spill=1)
#   dg$edg[dg$edg>1] <- 1
#   nrp <- nearest.raster.point(x=ppt$x,y=ppt$y,w=Wm)
#   ep.g.i <- dg$edg[cbind(nrp$row,nrp$col)]
#   
#   dg.dens <- spattemp.density(ppt,g,gam,tlim,TRUE,sres,tres,NA)
#   ts <- dg.dens$TT
#   qtlam <- getQt(ts,lam,tlim=tlim)
#   ep.gam.i <- getQt(tt,gam,tlim=tlim)
#   evalyx <- as.matrix(expand.grid(dg.dens$X,dg.dens$Y))
#   inn <- inside.owin(x=evalyx[,2],y=evalyx[,1],w=W)
#   
#   evalyx.redu <- evalyx[inn,] 
#   GN <- nrow(evalyx.redu)
#   dexi <- matrix(rep(dh$edg[which(inn)],tres),GN,tres)
#   deti <- matrix(rep(qtlam,GN),GN,tres,byrow=TRUE)
#   
#   
#   bs2temp <- bs1 <- bs2 <- matrix(NA,n,GN)
#   for(xi in 1:GN){
#     ev <- cbind(rep(evalyx.redu[xi,2],n),rep(evalyx.redu[xi,1],n))-ppdat
#     bs2temp[,xi] <- ep.g.i^(-1)*dmvnorm(ev,mean=c(0,0),sigma=diag(2)*(h^2+g^2))
#     bs2[,xi] <- ep.g.i^(-1)*dmvnorm(ev,mean=c(0,0),sigma=diag(2)*g^2)
#     bs1[,xi] <- ep.g.i^(-1)*dmvnorm(ev,mean=c(0,0),sigma=diag(2)*(0.5*h^2+g^2))
#   }
#   
#   bt2temp <- bt1 <- bt2 <- matrix(NA,n,tres)
#   b1c <- b2c <- matrix(NA,GN,tres)
#   for(tti in 1:tres){
#     et <- rep(ts[tti],n)-tt
#     bt2temp[,tti] <- ep.gam.i^(-1)*dnorm(et,0,sqrt(lam^2+gam^2))
#     bt2[,tti] <- ep.gam.i^(-1)*dnorm(et,0,gam)
#     bt1[,tti] <-  ep.gam.i^(-1)*dnorm(et,0,sqrt(0.5*lam^2+gam^2))
#     
#     b2temp.x.tti <- colSums(bs2temp*matrix(rep(bt2temp[,tti],GN),n,GN))
#     b2c[,tti] <- b2temp.x.tti*colSums(bs2*matrix(rep(bt2[,tti],GN),n,GN))
#     b1c[,tti] <- (1/(8*pi^1.5*h^2*lam))*colSums(bs1*matrix(rep(bt1[,tti],GN),n,GN)) + ((n-1)/n)*b2temp.x.tti^2
#   }		
#   
#   cella <- (dg.dens$X[2]-dg.dens$X[1])*(dg.dens$Y[2]-dg.dens$Y[1])*(dg.dens$TT[2]-dg.dens$TT[1])
#   
#   boot1 <- deti^(-2)*dexi^(-2)*b1c
#   boot2 <- -2*deti^(-1)*dexi^(-1)*b2c
#   boot3 <- sum(dg.dens$Z^2*cella,na.rm=TRUE)
#   
#   return(sum(n^(-2)*(boot1+boot2)*cella,na.rm=TRUE)+boot3)
# }
