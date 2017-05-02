tol.asy.ada <- function(f,g,beta,verbose=FALSE){
  
  if(verbose) cat("Initialising window...")
  if(!compatible(f$z,g$z)) stop("incompatible images 'f' and 'g' -- kernel estimates must be evaluated on identical domains")
  
  M <- Window(f$z)
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
  
  symmetric <- all(f$him==g$him,na.rm=TRUE)
  fq <- ifelse(is.null(f$q)||is.vector(f$q),NA,f$q)
  gq <- ifelse(is.null(g$q)||is.vector(g$q),NA,g$q)
  fqn <- is.na(fq)
  gqn <- is.na(gq)
  needq <- FALSE
  if(fqn||gqn){
    if(!symmetric){
      needq <- TRUE
    } else {
      if(fqn&&!gqn){
        fq <- gq
      } else if(gqn&&!fqn){
        gq <- fq
      } else {
        needq <- TRUE
      }
    }
  }
  
  qb <- function(oo,hfac,callno){
    hfp <- hfac*oo$him
    hfp$v[is.na(hfp$v)] <- max(hfp,na.rm=TRUE)
    hypoQ <- unique(quantile(hfp[inside],seq(0,1,beta),na.rm=TRUE))
    hypocen <- hypoQ[-length(hypoQ)]+diff(hypoQ/2)
    corrQ <- as.numeric(cut(as.vector(as.matrix(hfp)),breaks=hypoQ,include.lowest=TRUE))
    if(verbose){
      if(callno==1) cat("Convolving bandwidth-categorised kernels with window:\n --pass 1...")
      else cat(paste(" --pass ",callno,"...",sep=""))
    }
    qhz <- rep(NA,pres^2)
    for(i in 1:length(hypocen)){
      densX.ker <- dnorm(xcol.ker,sd=hypoQ[i])
      densY.ker <- dnorm(yrow.ker,sd=hypoQ[i])
      Kern <- outer(densY.ker,densX.ker,"*")*kerpixarea
      con <- fft(fM*fft(Kern), inverse=TRUE)/len.pad
      edg <- im(Mod(con[1:pres,1:pres]),xcol=oo$z$xcol,yrow=oo$z$yrow)
      qhz[which(corrQ==i)] <- as.vector(as.matrix(edg))[which(corrQ==i)]
    }
    if(verbose) cat("Done.\n")
    qhz[qhz>1] <- 1
    qhz <- im(matrix(qhz,pres,pres),xcol=oo$z$xcol,yrow=oo$z$yrow)
    qhz[!inside] <- NA
    return(qhz)
  }
  
  qd <- function(oo,hfac,callno){
    evalxy <- as.matrix(expand.grid(oo$z$xcol,oo$z$yrow))
    notin <- !inside.owin(x=evalxy[,1],y=evalxy[,2],w=M)
    h.hypo.mat <- as.matrix(hfac*oo$him)
    if(verbose){
      if(callno==1) cat("Processing pixel-by-pixel factors:\n --pass 1...")
      else cat(paste(" --pass ",callno,"...",sep=""))
    }
    qhz <- rep(NA,pres^2)
    for(i in 1:nrow(evalxy)){
      ht <- h.hypo.mat[which(oo$z$yrow==evalxy[i,2]),which(oo$z$xcol==evalxy[i,1])]
      if(is.na(ht)) next
      
      gxy <- ht^(-2)*(exp(-0.5*rowSums((cbind(evalxy[,1]-evalxy[i,1],evalxy[,2]-evalxy[i,2])/ht)^2))/(2*pi))
      gxy[notin] <- NA
      qhz[i] <- integral(im(matrix(gxy,pres,pres,byrow=TRUE),xcol=oo$z$xcol,yrow=oo$z$yrow))
    }
    if(verbose) cat("Done.\n")
    return(im(matrix(qhz,pres,pres,byrow=TRUE),xcol=oo$z$xcol,yrow=oo$z$yrow))
  }
  
  qob <- function(oo){
    if(!is.na(beta)){
      q <- oo$q
      q2 <- qb(oo,sqrt(0.5),callno)
      if(needq) q <- qb(oo,1,callno+1)
    } else {
      q <- oo$q
      q2 <- qd(oo,sqrt(0.5),callno)
      if(needq) q <- qd(oo,1,callno+1)
    }
    return(list(q=q,q2=q2))
  }
  
  if(verbose) cat("Done.\n")
  
  callno <- 1
  qqf <- qqg <- qob(f)
  if(!symmetric){
    callno <- 2 + needq
    qqg <- qob(g)
  }
  
  fk2 <- (1/(4*pi))*as.matrix(qqf$q2)
  gk2 <- (1/(4*pi))*as.matrix(qqg$q2)
  
  S1rzK <- (1/as.matrix(qqf$q)^2)*(2*fk2 + 0.5*fk2)
  S2rzK <- (1/as.matrix(qqg$q)^2)*(2*gk2 + 0.5*gk2)
  
  numerator <- as.matrix(log(f$z)-log(g$z))
  denominator <- sqrt(((S1rzK*f$gamma^2)/(length(f$h)*f$h0^2))+((S2rzK*g$gamma^2)/(length(g$h)*g$h0^2)))
  zstandard <- numerator/denominator
  P <- pnorm(zstandard,lower.tail=FALSE)
  return(list(num=numerator,den=denominator,zst=zstandard,p=im(P,xcol=f$z$xcol,yrow=f$z$yrow),s1=S1rzK,s2=S2rzK))
}