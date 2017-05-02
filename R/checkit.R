checkit <- function(h,str){
  if(!is.numeric(h)) stop(paste(str,"must be numeric"))
  if(length(h)>1){
    warning(paste(str,"must have length 1; using first element only"))
    h <- h[1]
  }
  if(h<=0) stop(paste(str,"must be positive"))
  return(h)
}