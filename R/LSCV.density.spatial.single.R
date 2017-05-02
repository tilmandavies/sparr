LSCV.density.spatial.single <- function(h,pp,res,edge){
  if(h<=0) return(NA)
  temp.dens <- density(pp,h,edge=edge,dimyx=res,positive=TRUE,diggle=FALSE)
  temp.int <- integral(temp.dens)
  temp.dens.pts <- density(pp,sigma=h,edge=edge,dimyx=res,at="points",positive=FALSE,leaveoneout=TRUE,diggle=FALSE)/temp.int
  temp.dens <- temp.dens/temp.int
  if(any(temp.dens.pts<=0)) return(NA) ## tiny bandwidth protector
  return(integral(temp.dens^2)-2*mean(temp.dens.pts))
}
