BAM.single <- function(h,edge,BP){
  cases <- BP$cas
  controls <- BP$con
  erode <- BP$erode
  xyin <- BP$xyin
  res <- BP$res
  kpa <- BP$kpa
  xk <- BP$xk
  yk <- BP$yk
  lp <- BP$lp
  fM <- BP$fM
  WM <- BP$WM
  fdd <- im(matrix(BP$fdd,res,res),xcol=WM$xcol,yrow=WM$yrow)
  gdd <- im(matrix(BP$gdd,res,res),xcol=WM$xcol,yrow=WM$yrow)
  
  qb <- function(hfac){
    hfp <- hfac*h
    densX.ker <- dnorm(xk,sd=h)
    densY.ker <- dnorm(yk,sd=h)
    Kern <- outer(densY.ker,densX.ker,"*")*kpa
    con <- fft(fM*fft(Kern), inverse=TRUE)/lp
    qhz <- im(Mod(con[1:res,1:res]),xcol=WM$xcol,yrow=WM$yrow)
    qhz[qhz>1] <- 1
    qhz[!WM$m] <- NA
    return(qhz)
  }
  
  qq <- qb(1)
  fd <- density.ppp(cases,sigma=h,dimyx=res,positive=TRUE,edge=FALSE)
  gd <- density.ppp(controls,sigma=h,dimyx=res,positive=TRUE,edge=FALSE)
  
  if(edge){
    fd <- fd/qq
    gd <- gd/qq
    qq <- as.vector(as.matrix(qq))
    qq2 <- as.vector(as.matrix((1/(4*pi))*qb(sqrt(0.5))))
    qq[!xyin] <- NA
    qq2[!xyin] <- NA
    qqi <- im(matrix(qq,res,res),xcol=WM$xcol,yrow=WM$yrow)
    qq2i <- im(matrix(qq2,res,res),xcol=WM$xcol,yrow=WM$yrow)
    rk <- integral(qq2i)/(qqi^2*h^2) ### now using 'image' qq
  } else {
    rk <- rep(1/(4*pi),res^2)
    rk[!xyin] <- NA
    rk <- im(matrix(rk,res,res),xcol=WM$xcol,yrow=WM$yrow)
  }
  
  fd <- fd/integral(fd)
  gd <- gd/integral(gd)
  
  # RK <- function(xs,ys,h,W,ca,xyin){
  #   qres <- qhzfunc(xs,ys,"gaus",W,h)
  #   qres$qhz[!xyin] <- NA
  #   qres$qhz_sq[!xyin] <- NA
  #   return(sum(qres$qhz_sq*ca,na.rm=T)/(qres$qhz^2*h^2))
  # }
  
  return(zeta1(fd,gd,npoints(cases),npoints(controls),rk)/h^2 + 0.5*h^4*zeta2(fd,gd,fdd,gdd))
  
}


zeta1 <- function(fh,gh,n1,n2,rk){
  return(integral(rk/fh)/n1+integral(rk/gh)/n2)
}

zeta2 <- function(fh,gh,fdd,gdd){
  return(integral(fdd^2/fh^2)/2-integral((fdd*gdd)/(fh*gh))+integral(gdd^2/gh^2)/2)
}
