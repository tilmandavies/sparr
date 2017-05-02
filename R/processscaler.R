processscaler <- function(s,pp){
  if(is.numeric(s)){
    if(length(s)>1){
      s <- s[1]
      warning("'scaler' if numeric must be of length 1. Using first value only.")
    }
    if(s<=0) stop("'scaler' must be positive")
  } else if(is.character(s)){
    if(length(s)>1){
      s <- s[1]
    }
    s <- switch(s,IQR=mean(c(IQR(pp$x)/1.34,IQR(pp$y)/1.34)),
                sd=mean(c(sd(pp$x),sd(pp$y))),
                var=sqrt(mean(c(var(pp$x),var(pp$y)))),
                silverman=min(mean(c(sd(pp$x),sd(pp$y))),mean(c(IQR(pp$x)/1.34,IQR(pp$y)/1.34))),NA)
    if(is.na(s)){
      stop("'scaler' character string must be one of \"silverman\", \"IQR\", \"sd\", or \"var\"")
    }
  } else {
    stop("Invalid 'scaler' type")
  }
  return(s)
}
