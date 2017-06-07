LIK.density.spatial.single <- function(h,pp,res,edge){
  if(h<=0) return(NA)
  temp.dens.pts <- density(pp,sigma=h,edge=edge,dimyx=res,at="points",positive=FALSE,leaveoneout=TRUE,diggle=FALSE)/npoints(pp)
  if(any(temp.dens.pts<=0)) return(log(min(temp.dens.pts))) ## tiny bandwidth protector
  return(mean(log(temp.dens.pts)))
}
