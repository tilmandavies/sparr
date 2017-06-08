processnstar <- function(n,pp){
  if(is.numeric(n)){
    if(length(n)>1){
      n <- n[1]
      warning("'nstar' if numeric must be of length 1. Using first value only.")
    }
    if(n<=0) stop("'nstar' must be positive")
  } else if(is.character(n)){
    if(length(n)>1){
      n <- n[1]
    }
    if(n=="npoints"){
      n <- npoints(pp)
    } else if(n=="geometric"){
      pm <- marks(pp)
      if(is.null(pm)||(!is.factor(pm))){
        n <- npoints(pp)
        warning("Using 'nstar'=\"geometric\" requires a 'ppp' object with factor-valued marks. Changing to \"npoints\"")
      } else {
        n <- prod(as.numeric(table(pm)))^(1/length(levels(pm)))
      }
    } else {
      stop("'nstar' character string must be one of \"npoints\" or \"geometric\"")
    }
  } else {
    stop("Invalid 'nstar' type")
  }
  return(n)
}

processnstar.st <- function(n,pp){
  if(is.numeric(n)){
    if(length(n)>1){
      n <- n[1]
      warning("'nstar' if numeric must be of length 1. Using first value only.")
    }
    if(n<=0) stop("'nstar' must be positive")
  } else if(is.character(n)){
    if(length(n)>1){
      n <- n[1]
    }
    if(n=="npoints"){
      n <- npoints(pp)
    } else if(n=="geometric"){
      n <- npoints(pp)
      warning("'nstar' character string cannot be \"geometric\" for spatiotemporal estimates -- using \"npoints\"")
    } else {
      stop("'nstar' character string only permitted to be \"npoints\"")
    }
  } else {
    stop("Invalid 'nstar' type")
  }
  return(n)
}
