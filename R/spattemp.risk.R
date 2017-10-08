#' @export
spattemp.risk <- function(f,g,log=TRUE,tolerate=FALSE,finiteness=TRUE,verbose=TRUE){
  if(!inherits(f,"stden")) stop("'f' must be of class 'stden' arising from a call to 'spattemp.density'")
  
  gst <- inherits(g,"stden")
  gbi <- inherits(g,"bivden")
  if(!(gst||gbi)) stop("'g' must be of class 'stden' or 'bivden'")
  
  fse <- all(f$qs==1)
  fte <- all(f$qt==1)
  
  if(verbose) message("Calculating ratio...", appendLF=FALSE)
  fres <- f$z
  flen <- length(fres)
  rr <- rrc <- list()
  if(gst){
    gse <- all(g$qs==1)
    gte <- all(g$qt==1)
    if((fse!=gse)||(fte!=gte)) stop("edge-correction for 'f' and 'g' must be consistent")
    
    gres <- g$z
    if(flen!=length(gres)) stop("incompatible temporal domains... 'f' and 'g' must be evaluated at identical timestamps")
    fn <- as.numeric(names(fres))
    gn <- as.numeric(names(gres))
    if(any(fn!=gn)) stop("incompatible temporal domains... 'f' and 'g' must be evaluated at identical timestamps")
    if(!compatible(fres[[1]],gres[[1]])) stop("incompatible images in 'f' and 'g'... kernel estimates must be evaluated on identical spatial domains")

    # if(positive){
    #   fres <- lapply(fres,posifybivden)
    #   gres <- lapply(gres,posifybivden)
    #   f$z.cond <- lapply(f$z.cond,posifybivden)
    #   g$z.cond <- lapply(g$z.cond,posifybivden)
    # }
    
    for(i in 1:flen){
      rr[[i]] <- suppressWarnings(log(fres[[i]])-log(gres[[i]]))
      rrc[[i]] <- suppressWarnings(log(f$z.cond[[i]])-log(g$z.cond[[i]]))
    }
    

    if(verbose) message("Done.")
    if(tolerate){
      if(verbose) message("Calculating pooled estimate for tolerance...", appendLF=FALSE)
      fs <- f
      gs <- g
      marks(fs$pp) <- NULL
      marks(gs$pp) <- NULL
      pdat <- suppressWarnings(superimpose(fs$pp,gs$pp))
      marks(pdat) <- c(marks(f$pp),marks(g$pp))
      hpool <- sqrt(prod(c(f$h,g$h)))
      lpool <- sqrt(prod(c(f$lambda,g$lambda)))
      pooled <- spattemp.density(pdat,h=hpool,lambda=lpool,tlim=f$tlim,sedge=ifelse(fse,"none","uniform"),tedge=ifelse(fte,"none","uniform"),sres=nrow(f$spatial.z),tres=flen,verbose=FALSE)
      if(verbose) message("Done.")
    }

  } else {
    if(!compatible(fres[[1]],g$z)) stop("incompatible images in 'f' and 'g'... kernel estimates must be evaluated on identical spatial domains")
    
    gse <- is.null(g$q)
    if(fse!=gse) stop("edge-correction for 'f' and 'g' must be consistent")
    
    g$z <- g$z/integral(g$z)
    
    # if(positive){
    #   fres <- lapply(fres,posifybivden)
    #   g$z <- posifybivden(g$z)
    #   f$z.cond <- lapply(f$z.cond,posifybivden)
    # }
    
    for(i in 1:flen){
      rr[[i]] <- suppressWarnings(log(fres[[i]])-log(g$z)+log(diff(f$tlim)))
      rrc[[i]] <- suppressWarnings(log(f$z.cond[[i]])-log(g$z))
    }
    
    if(verbose) message("Done.")
    pooled <- NULL
  }
  
  if(finiteness&&log){
    if(verbose) message("Ensuring finiteness...\n   --joint--")
    rr <- lapply(rr,fbound)
    if(verbose) message("   --conditional--")
    rrc <- lapply(rrc,fbound)
    if(verbose) message("Done.")
  }
  
  
  ps <- psc <- NULL
  if(tolerate){
    if(verbose) message("Calculating tolerance contours...")
    vars <- tol.asy.st(f,g,pooled,verbose)
    ps <- psc <- list()
    for(i in 1:flen){
      Z <- as.matrix(rr[[i]])/sqrt(vars$v[[i]])
      ps[[i]] <- im(pnorm(Z,lower.tail=FALSE),xcol=f$spatial.z$xcol,yrow=f$spatial.z$yrow)
      Zc <- as.matrix(rrc[[i]])/sqrt(vars$vc[[i]])
      psc[[i]] <- im(pnorm(Zc,lower.tail=FALSE),xcol=f$spatial.z$xcol,yrow=f$spatial.z$yrow)
    }
    if(verbose) message("Done.")
  }
  
  if(!log){
    for(i in 1:flen){
      rr[[i]] <- exp(rr[[i]])
      rrc[[i]] <- exp(rrc[[i]])
    }
  }
  
  names(rr) <- names(rrc) <- names(f$z)
  if(!is.null(ps)) names(ps) <- names(psc) <- names(f$z)
  result <- list(rr=rr,rr.cond=rrc,P=ps,P.cond=psc,f=f,g=g,tlim=f$tlim)
  class(result) <- "rrst"
  return(result)
}

fbound <- function(x){
  fnt <- is.finite(as.matrix(x))
  if(!any(fnt)){
    x[] <- 0
  } else {
    imr <- range(x[fnt])
    x[x==Inf] <- imr[2]
    x[x==-Inf] <- imr[1]
  }
  return(x)
}