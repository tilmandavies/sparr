KO <- function(xy){
  res <- exp(-0.5*rowSums(xy^2))/(2*pi)
  return((2*res-xy[,1]^2*res-xy[,2]^2*res)^2)
}