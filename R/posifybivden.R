posifybivden <- function(x, eps=.Machine$double.xmin) {
  force(eps)
  if(is.im(x)) return(eval.im(pmax(eps, x)))
  if(is.numeric(x)) return(pmax(eps, x))
  return(x)
}
