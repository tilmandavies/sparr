gethcats <- function(h,breaks=NULL,step=0.05){
  if(is.null(breaks)) breaks <- unique(as.numeric(quantile(h,seq(0,1,step))))
  hc <- as.numeric(cut(h,breaks=breaks,include.lowest=TRUE,right=TRUE))
  hcen <- breaks[-length(breaks)]+diff(breaks)/2
  return(hcen[hc])
}
