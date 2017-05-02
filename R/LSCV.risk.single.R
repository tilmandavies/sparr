LSCV.risk.single <- function(h,cases,controls,res,edge,hazey){
  if(h<=0) return(NA)
  
  temp.case <- density.ppp(cases,sigma=h,edge=edge,dimyx=res,positive=TRUE)
  temp.con <- density.ppp(controls,sigma=h,edge=edge,dimyx=res,positive=TRUE)
  tcase.int <- integral(temp.case)
  tcon.int <- integral(temp.con)
  temp.case <- temp.case/tcase.int
  temp.con <- temp.con/tcon.int
  
  if(any(is.infinite(as.matrix(temp.case/temp.con)))) return(NA)  ## pre-fail for infinite rr cells - both HAZE and KELDIG
  
  temp.case.pts <- density.ppp(cases,sigma=h,edge=edge,dimyx=res,at="points",leaveoneout=TRUE,positive=TRUE)/tcase.int
  temp.con.pts <- density.ppp(controls,sigma=h,edge=edge,dimyx=res,at="points",leaveoneout=TRUE,positive=TRUE)/tcon.int
  
  caseatcon <- safelookup(temp.case,controls,warn=FALSE)
  conatcase <- safelookup(temp.con,cases,warn=FALSE)

  if(any(temp.case.pts<=0)||any(temp.con.pts<=0)||any(caseatcon<=0)||any(conatcase<=0)) return(NA) ## tiny bandwidth protector
  
  if(!hazey) result <- 2*mean(log(caseatcon/temp.con.pts)/temp.con.pts) - 2*mean(log(temp.case.pts/conatcase)/temp.case.pts) - integral((log(temp.case)-log(temp.con))^2)
  else result <- mean((caseatcon/temp.con.pts)^2)-2*mean(temp.case.pts/conatcase)
  
  return(result)
}
