available.h0 <- function(...){
  unpacked <- list(...)
  
  cls <- lapply(unpacked,function(x) inherits(x,"msden"))
  if(!all(unlist(cls))) stop("function arguments must all be of class \"msden\", arising from a call to 'multiscale.density'")
  
  lo <- sapply(unpacked,function(x) x$h0range[1])
  up <- sapply(unpacked,function(x) x$h0range[2])
  rng <- c(max(lo),min(up))
  
  if(rng[1]>=rng[2]) stop("incompatible 'h0range' components -- check bandwidth scales")
  return(rng)
}