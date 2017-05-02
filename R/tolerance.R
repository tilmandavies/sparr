tolerance <- function(rs, method = c("ASY", "MC"), ref.density = NULL, beta = 0.025, ITER = 100, parallelise = NULL, verbose = TRUE, ...){
  if(!inherits(rs,"rrs")) stop("'rs' argument must be of class \"rrs\"")
  
  meth <- method[1]
  adaf <- !is.null(rs$f$him)
  adag <- !is.null(rs$g$him)
  ada <- adaf + adag
  if(ada==1) stop("case/control smoothing regimens (fixed or adaptive) must be identical")
  
  if(meth=="ASY"){
    if(ada==2){
      psurf <- tol.asy.ada(rs$f,rs$g,beta,verbose)$p
    } else {
      if(is.null(ref.density)) stop("must supply 'ref.density' for fixed-bandwidth asymptotic tolerance contours")
      if((!inherits(ref.density,"bivden"))&&(!inherits(ref.density,"im"))) stop("'ref.density' must be of class \"bivden\" or \"im\"")
      if(is.im(ref.density)){
        ref.density <- list(z=ref.density,q=NULL)
      #'  was:
      #   rq <- ref.density
      # 	rq$v[!is.na(rq$v)] <- 1
      # 	ref.density <- list(z=ref.density,q=rq)
      }
      ref.density$z <- ref.density$z/integral(ref.density$z)
      if(!compatible(rs$f$z,rs$g$z,ref.density$z)) stop("incompatible 'ref.density'... must be evaluated on domain identical to case/control densities")
      psurf <- tol.asy.fix(rs$f,rs$g,ref.density,verbose)$p
    }
  } else if(meth=="MC"){
    ITER <- checkit(ITER,"'ITER'")
    if(ada==2){
      psurf <- im(tol.mc.ada(rs,round(ITER),parallelise,verbose,davies.baddeley=beta,...),xcol=rs$rr$xcol,yrow=rs$rr$yrow)
    } else {
      psurf <- im(tol.mc.fix(rs,round(ITER),parallelise,verbose,...),xcol=rs$rr$xcol,yrow=rs$rr$yrow)
    }
  } else stop("invalid 'method' argument")
  
  return(psurf)
}
  