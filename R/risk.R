risk <- function(f, g = NULL, log = TRUE, h0 = NULL, hp = h0, adapt = FALSE,
                  tolerate = FALSE, doplot = FALSE,
                  pilot.symmetry = c("none","f","g","pooled"), epsilon = 0,
                  verbose = TRUE, ...){

  if(is.null(g)){
    if(!inherits(f,"ppp")) stop("'f' must be an object of class 'ppp' if 'g' unsupplied")
    fm <- marks(f)
    if(!is.factor(fm)) marks(f) <- fm <- factor(fm)
    if(nlevels(fm)!=2) stop("'f' marks must be dichotomous if 'g' unsupplied")
    fs <- split(f)
    f <- fs[[1]]
    g <- fs[[2]]
  } else {
    fc <- class(f)
    gc <- class(g)
    if(!all(fc==gc)) stop("'f' and 'g' must be of identical class")
    if(!(inherits(f,"ppp")||inherits(f,"bivden"))) stop("'f' and 'g' must be of class 'ppp' or 'bivden'")
  }
  
  epsi <- epsilon[1]
  if(epsi<0) stop("invalid 'epsilon'; must be scalar and non-negative")
  
  if(inherits(f,"ppp")){
    if(!identical_windows(Window(f),Window(g))) stop("study windows for 'f' and 'g' must be identical")
    
    pooled <- suppressWarnings(superimpose(f,g))
    if(is.null(h0)) h0 <- OS(pooled)
    
    if(length(h0)==1){
      h0f <- h0g <- checkit(h0[1],"'h0[1]'")
    } else {
      h0f <- checkit(h0[1],"'h0[1]'")
      h0g <- checkit(h0[2],"'h0[2]'")
    }
    
    if(!adapt){
      if(verbose) cat("Estimating case and control densities...")
      fd <- bivariate.density(f,h0=h0f,adapt=FALSE,...)
      gd <- bivariate.density(g,h0=h0g,adapt=FALSE,...)
      if(verbose) cat("Done.\n")
    } else {
      
      if(is.null(hp)) hp <- c(h0f,h0g)
      
      if(length(hp)==1){
        hfp <- hgp <- checkit(hp[1],"'hp[1]'")
      } else {
        hfp <- checkit(hp[1],"'hp[1]'")
        hgp <- checkit(hp[2],"'hp[2]'")
      }
      
      pilot.symmetry <- pilot.symmetry[1]
      pilotdata <- switch(pilot.symmetry,none=1,f=f,g=g,pooled=pooled,NA)
      if(any(is.na(pilotdata))) stop("invalid 'pilot.symmetry' argument")
      
      if(verbose) cat("Estimating pilot(s)...")
      if(pilot.symmetry=="none"){
        fp <- bivariate.density(f,h0=hfp,adapt=FALSE,...)
        gp <- bivariate.density(g,h0=hgp,adapt=FALSE,...)
      } else {
        fp <- gp <- bivariate.density(pilotdata,h0=hfp[1],adapt=FALSE,...)
      }
      
      fgeo <- log(posifybivden(safelookup(fp$z,f,warn=FALSE))^(-0.5))
      ggeo <- log(posifybivden(safelookup(gp$z,g,warn=FALSE))^(-0.5))
      gam <- exp(npoints(pooled)^(-1)*(sum(fgeo)+sum(ggeo)))
      
      if(verbose) cat("Done.\n")
      
      if(verbose) cat("Estimating case density...")
      fd <- bivariate.density(f,h0=h0f,adapt=TRUE,pilot.density=fp$z,gamma.scale=gam,verbose=FALSE,...)
      if(verbose) cat("Done.\nEstimating control density...")
      gd <- bivariate.density(g,h0=h0g,adapt=TRUE,pilot.density=gp$z,gamma.scale=gam,verbose=FALSE,...)
      if(verbose) cat("Done.\n")
    }
  } else {
    if(!compatible(f$z,g$z)) stop("incompatible images in 'f' and 'g'... kernel estimates must be evaluated on identical domains")
    fd <- f
    gd <- g
    fda <- is.na(fd$gamma)||is.na(fd$geometric)
    gda <- is.na(gd$gamma)||is.na(gd$geometric)
    adapt <- switch(as.character(fda+gda),"0"=TRUE,"2"=FALSE,NA)
    if(is.na(adapt)) stop("'f' and 'g' smoothed differently... must both be either fixed or adaptive")
  }
  
  eg <- epsi*max(gd$z)
  #rr <- (fd$z+eg)/(gd$z+eg)
  #if(log) rr <- log(rr)
  
  if(log) suppressWarnings(rr <- log(fd$z+eg) - log(gd$z+eg))
  else rr <- (fd$z+eg)/(gd$z+eg)
  
  ps <- NULL
  if(tolerate){
    if(verbose) cat("Calculating tolerance contours...")
    if(adapt) ps <- tol.asy.ada(fd,gd,0.025,verbose=FALSE)$p
    else ps <- tol.asy.fix(fd,gd,gd,verbose=FALSE)$p
    if(verbose) cat("Done.\n")
  }
  
  if(doplot){
    plot.im(rr,main="",box=FALSE,ribargs=list(box=TRUE))
    axis(1)
    axis(2)
    box(bty="l")
    plot(Window(fd$pp),add=TRUE)
    if(!is.null(ps)) contour(fd$z$xcol,fd$z$yrow,t(as.matrix(ps)),levels=0.05,add=TRUE)
  }
  
  result <- list(rr=rr,f=fd,g=gd,P=ps)
  class(result) <- "rrs"
  
  return(result)
}
