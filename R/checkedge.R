checkedge <- function(edge,v=1){
  if(!is.vector(edge)) stop("'edge' argument must be a vector")
  if(!is.character(edge)) stop("'edge' argument must be a character string")
  etype <- edge[1]
  if(v==1){
    if(!any(etype==c("uniform","diggle","none"))) stop("invalid 'edge' argument")
  } else {
    if(!any(etype==c("uniform","none"))) stop("invalid 'edge' argument")
  }
  return(etype)
}
