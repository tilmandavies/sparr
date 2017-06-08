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


processscaler.st <- function(s,pp,tt){
  if(is.numeric(s)){
    if(length(s)==1){
      s <- c(s,s)
    } else if(length(s)>1){
      s <- s[1:2] 
    } else {
      stop("If numeric, 'scaler' must be of length 1 or 2")
    }
    if(s<=0) stop("'scaler' must be wholly positive")
  } else if(is.character(s)){
    if(length(s)>1){
      s <- s[1]
    }
    s <- switch(s,IQR=c(mean(c(IQR(pp$x)/1.34,IQR(pp$y)/1.34)),IQR(tt)/1.34),
                sd=c(mean(c(sd(pp$x),sd(pp$y))),sd(tt)),
                var=c(sqrt(mean(c(var(pp$x),var(pp$y)))),sd(tt)),
                silverman=c(min(mean(c(sd(pp$x),sd(pp$y))),mean(c(IQR(pp$x)/1.34,IQR(pp$y)/1.34))),min(sd(tt),IQR(tt)/1.34)),NA)
    if(any(is.na(s))){
      stop("'scaler' character string must be one of \"silverman\", \"IQR\", \"sd\", or \"var\"")
    }
  } else {
    stop("Invalid 'scaler' type")
  }
  return(s)
}
