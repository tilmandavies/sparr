#' @export
BOOT.density <- function(pp,hlim=NULL,eta=NULL,type=c("fixed","adaptive"),hp=NULL,
                         edge=c("uniform","none"),ref.density=NULL,resolution=64,
                         rmdiag=TRUE,sim.adapt=list(N=50,B=100,dimz=64,objective=FALSE),
                         parallelise=NA,verbose=TRUE,...){
  
  if(!inherits(pp,"ppp")) stop("'pp' must be of spatstat class \"ppp\"; see ?ppp")
  
  n <- npoints(pp)
  W <- Window(pp)
  WM <- as.mask(W,dimyx=rep(resolution,2))
  
  # set h-limits if unsupplied
  if(is.null(hlim)){
    ppu <- pp
    marks(ppu) <- NULL
    md <- min(nndist(unique(ppu)))
    hlim <- c(md,max(md*50,min(diff(W$xrange),diff(W$yrange))/6))
  } else {
    hlim <- checkran(hlim,"'hlim'")
  }
  
  edg <- checkedge(edge,v=0)
  
  if(!is.na(parallelise)){
    if(!is.numeric(parallelise)) stop("'parallelise' must be numeric")
    if(is.null(parallelise)) parallelise <- NA
    parallelise <- round(parallelise[1])
  }
  
  typ <- type[1]
  if(typ=="fixed"){
    
    if(verbose) cat("Initialising...")
    
    # unleash the grid within
    inside <- WM$m
    inn <- which(as.vector(inside))
    evalyx <- as.matrix(expand.grid(WM$yrow,WM$xcol))
    evalyx.redu <- evalyx[inn,]
    GN <- length(inn)
    
    edg <- edg=="uniform"
    
    # get reference density and associated diggle edge factors
    if(is.null(ref.density)){
      
      if(is.null(eta)){
        eta <- OS(pp)
      } else {
        eta <- checkit(eta,"'eta'")
      }
      
      d.eta <- density(pp,eta,edge=edg,positive=TRUE,dimyx=resolution,spill=1)
      if(edg){
        d.eta$edg[d.eta$edg>1] <- 1
        d.etadens <- d.eta$raw/d.eta$edg
        d.etaint <- integral(d.etadens)
        d.etadens <- as.matrix(d.etadens)[inn]/d.etaint
        epsilon.eta <- safelookup(d.eta$edg,pp,warn=FALSE)
      } else {
        d.etadens <- as.matrix(d.eta$raw/integral(d.eta$raw))[inn]
        epsilon.eta <- rep(1,n)
      }
      
    } else {
      if(!inherits(ref.density,"bivden")) stop("'ref.density' must be of class \"bivden\"; see ?bivden")
      
      eta <- ref.density$h0
      if(!compatible(as.im(WM),ref.density$z)) stop("'ref.density' must be evaluated on identical spatial domain as 'Window(pp)' given 'resolution'")
      if(!all(ref.density$h==eta)) stop("'ref.density' must be a fixed-bandwidth density estimate when type = \"fixed\"")
      if(edg&&is.null(ref.density$q)) stop("'ref.density' edge-correction must exist if edge = \"uniform\"")
      
      if(is.null(ref.density$q)){
        epsilon.eta <- rep(1,n)
      } else if(is.vector(ref.density$q)){
        epsilon.eta <- ref.density$q
      } else {
        epsilon.eta <- safelookup(ref.density$q,pp,warn=FALSE)
      }
      
      d.etadens <- ref.density$z/integral(ref.density$z)
      d.etadens <- as.matrix(d.etadens)[inn]
    }
    
    if(edg){
      # window dressing
      use_fftw <- fftw_available()
      ifft_scale <- WM$xstep*WM$ystep/(4*resolution^2)
      Mpad <- matrix(0,2*resolution,2*resolution)
      Mpad[1:resolution,1:resolution] <- inside
      fM <- fft2d(Mpad,fftw=use_fftw)
    } else {
      ifft_scale <- use_fftw <- fM <- NA
    }
    
    # fM <- fft2d(Mpad,fftw=use_fftw)
    # fK.eta <- kernel2d_fft(eta,WM$xstep,WM$ystep,resolution)
    # fK.con <- fft2d(fM*fK.eta,inverse=TRUE,use_fftw)[1:resolution,1:resolution]
    # edg.eta <- Mod(fK.con)*ifft_scale
    # edg.eta[edg.eta>1] <- 1  
    # edg.eta <- im(matrix(edge.eta,resolution,resolution),xcol=WM$xcol,yrow=WM$yrow)
    # epsilon.eta <- safelookup(edg.eta,pp,warn=FALSE)
  
    # remove expensive constant (just when rmdiag==TRUE) from 'boot.opt.spatial'; supply as argument
    if(!rmdiag) boot3 <- d.etadens^2
    else boot3 <- n^(-2)*sum(outer(1:n,1:n,function(i,j) epsilon.eta[i]^(-1)*epsilon.eta[j]^(-1)*kernel2d(pp$x[i]-pp$x[j],pp$y[i]-pp$y[j],sqrt(2*eta^2)))[-seq(1,n^2,n+1)])
  
    if(verbose) cat("Done.\nSearching for optimal h in ",prange(hlim),"...",sep="")
    result <- optimise(boot.opt.spatial.fix,interval=hlim,rmdiag=rmdiag,edg=edg,WM=WM,
                       resolution=resolution,fM=fM,ifft_scale=ifft_scale,
                       inn=inn,GN=GN,evalyx.redu=evalyx.redu,pp=pp,
                       epsilon.eta=epsilon.eta,eta=eta,nn=n,boot3=boot3,
                       use_fftw=use_fftw,parallelise=parallelise)$minimum
    if(verbose) cat("Done.\n")
  
  } else if(typ=="adaptive"){
    if(verbose) cat("Initialising...")
    
    if(is.null(ref.density)){
      if(is.null(eta)){
        eta <- OS(pp)
      } else {
        eta <- checkit(eta,"'eta'")
      }
      if(is.null(hp)){
        hp <- eta
      } else {
        hp <- checkit(hp,"'hp'")
      }
      d.eta <- density(pp,eta,edge=(edg=="uniform"),positive=TRUE,dimyx=resolution)
    } else {
      if(!inherits(ref.density,"bivden")) stop("'ref.density' must be of class \"bivden\"; see ?bivden")
      if(!compatible(as.im(WM),ref.density$z)) stop("'ref.density' must be evaluated on identical spatial domain as 'Window(pp)' given 'resolution'")
      if(is.null(hp)){
        hp <- ifelse(is.null(ref.density$hp),ref.density$h0,ref.density$hp)
      } else {
        hp <- checkit(hp,"'hp'")
      }
      d.eta <- ref.density$z
    }
    
    d.eta <- d.eta/integral(d.eta)
    
    if(is.null(sim.adapt)) sim.adapt <- list()
    if(is.null(sim.adapt$N)) sim.adapt$N <- 50
    if(is.null(sim.adapt$B)) sim.adapt$B <- 100
    if(is.null(sim.adapt$dimz)) sim.adapt$dimz <- 64
    if(is.null(sim.adapt$objective)) sim.adapt$objective <- FALSE
    
    hhash <- mean(hlim)
    h0fac <- hlim/hhash
  
    if(verbose) cat("Done.\nSearching for optimal h0 in ",prange(hlim),":\n",sep="")
    resultfull <- boot.opt.spatial.adapt(pp,h0ref=hhash,h0fac=h0fac,hp=hp,
                                         edg=edg,refden=d.eta,N=sim.adapt$N,
                                         B=sim.adapt$B,res=resolution,dimz=sim.adapt$dimz,
                                         verbose=verbose,parallelise=parallelise,...)
    
    
    result <- resultfull$h
    if(sim.adapt$objective) result <- resultfull$mat
    
  } else {
    stop("invalid 'type'")
  }
  
  return(result)
}









