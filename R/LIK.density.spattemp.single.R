LIK.density.spattemp.single <- function(bands,pp,tt,tlim,xyin,xys,sedge,tedge,parallelise,verbose){
  if(any(bands<=0)) return(NA)
  if(verbose) cat("h =",bands[1],"\b; lambda =",bands[2],"\n")
  h <- bands[1]
  lam <- bands[2]
  
  temp.dens.pts <- spattemp.LOO(pp,tt,h,lam,tlim,xyin,xys,sedge,tedge,parallelise=parallelise)
  if(any(temp.dens.pts<=0)) return(log(min(temp.dens.pts)))
  return(-mean(log(temp.dens.pts)))
}
