LSCV.density.spatial.single <- function(h,pp,res,edge){
  if(h<=0) return(NA)
  temp.dens <- density(pp,h,edge=edge,dimyx=res,positive=TRUE,diggle=FALSE)
  # if(any(temp.dens<=0)) return(NA)
  temp.int <- integral(temp.dens)
  temp.dens.pts <- density(pp,sigma=h,edge=edge,dimyx=res,at="points",positive=FALSE,leaveoneout=TRUE,diggle=FALSE)/temp.int
  temp.dens <- temp.dens/temp.int
  t2int <- integral(temp.dens^2)
  if(any(temp.dens.pts<=0)) return(t2int) ## tiny bandwidth protector
  return(t2int-2*mean(temp.dens.pts))
}
