LSCV.density.spatial.single <- function(h,pp,res,edge,za){
  if(h<=0) return(NA)
  temp.dens <- density(pp,h,edge=edge,dimyx=res,positive=FALSE,diggle=FALSE)
  temp.int <- integral(temp.dens)
  temp.dens.pts <- density(pp,sigma=h,edge=edge,dimyx=res,at="points",positive=FALSE,leaveoneout=TRUE,diggle=FALSE)/temp.int
  temp.dens <- temp.dens/temp.int
  t2int <- integral(temp.dens^2)
  
  ## tiny bandwidth protector action
  if(za==-1){
    if(any(temp.dens<=0)) return(Inf)
  }
  
  if(any(temp.dens.pts<=0)){ #was: return(t2int)
    if(za==2){
      temp.dens.pts[temp.dens.pts<=0] <- min(temp.dens.pts[temp.dens.pts>0])
    } else if(za==1){
      temp.dens.pts <- posifybivden(temp.dens.pts)
    } else {
      return(Inf)
    }
  }  
  return(t2int-2*mean(temp.dens.pts))
}  