subpix <- function(h,plocs,hc,WM,weights=NULL){
  id <- which(hc==h)
  res <- WM$dim[1]
  ploc.tempr <- plocs$row[id]
  ploc.tempc <- plocs$col[id]
  rf <- factor(ploc.tempr,levels=1:res)
  cf <- factor(ploc.tempc,levels=1:res)
  if(is.null(weights)) ta <- table(row=rf,col=cf)
  else ta <- tapply(weights[id],list(row=rf,col=cf),sum)
  ta[is.na(ta)] <- 0
  return(im(ta,xcol=WM$xcol,yrow=WM$yrow,xrange=WM$xrange,yrange=WM$yrange))
}