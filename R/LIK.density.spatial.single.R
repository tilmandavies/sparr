LIK.density.spatial.single <- function(h,pp,res,edge,za){
  if(h<=0) return(NA)
  temp.dens.pts <- density(pp,sigma=h,edge=edge,dimyx=res,at="points",positive=FALSE,leaveoneout=TRUE,diggle=FALSE)/npoints(pp)
  
  ## tiny bandwidth protector action
  if(za==-1){
    dtest <- density(pp,sigma=h,edge=edge,dimyx=res,diggle=FALSE)
    if(any(dtest<=0)) return(-Inf)
  }
  
  if(any(temp.dens.pts<=0)){ 
    if(za==2){
      temp.dens.pts[temp.dens.pts<=0] <- min(temp.dens.pts[temp.dens.pts>0])
    } else if(za==1){
      temp.dens.pts <- posifybivden(temp.dens.pts)
    } else {
      return(-Inf)
    }
  }
  return(mean(log(temp.dens.pts)))
}
