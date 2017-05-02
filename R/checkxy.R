checkxy <- function(xy){
  if(!is.list(xy)) stop("'xy' must be a list")
  if(is.null(xy$x)||is.null(xy$y)) stop("'xy' have vector members 'x' and 'y'")
  if(length(xy$x)!=length(xy$y)) stop("'xy' must be a list with numeric vectors 'x' and 'y' of equal length")
  if((!is.numeric(xy$x))||(!is.numeric(xy$y))) stop("'xy' must have numeric vector members 'x' and 'y'")
  return(xy)
}