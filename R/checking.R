checkdb <- function(db){
  if(!is.vector(db)) stop("'davies.baddeley' argument must be a vector")
  if(length(db)!=3&&length(db)!=1) stop("invalid 'davies.baddeley' argument length")
  if(!is.numeric(db)) stop("'davies.baddeley' argument must be numeric")
  db1 <- db[1]<=0||db[1]>=1
  db23 <- FALSE
  if(length(db)==3){
    db23 <- db[2]<=0||db[2]>=1||db[3]<=1
    dbres <- c(db[1:2],ceiling(db[3]))
  } else {
    dbres <- c(db,db,NA)
  }
  if(db1||db23) stop("one or more invalid 'davies.baddeley' values")
  return(dbres)
}

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

checkit <- function(h,str){
  if(!is.numeric(h)) stop(paste(str,"must be numeric"))
  if(length(h)>1){
    warning(paste(str,"must have length 1; using first element only"))
    h <- h[1]
  }
  if(h<=0) stop(paste(str,"must be positive"))
  return(h)
}

checkran <- function(ran,nm){
  if(!is.vector(ran)) stop(paste("'",nm,"' must be a vector",sep=""))
  if(!is.numeric(ran)) stop(paste("'",nm,"' must be numeric",sep=""))
  if(length(ran)!=2) stop(paste("'",nm,"' must have length 2",sep=""))
  if(ran[2]<=ran[1]) stop(paste("'",nm,"[1]' must be < '",nm,"[2]'",sep=""))
  if(any(ran<=0)) stop(paste("'",nm,"' must be wholly positive",sep=""))
  return(ran)
}

checkranin <- function(ran,vals,nm1){
  if(!is.vector(ran)) stop(paste("'",nm1,"' must be a vector",sep=""))
  if(!is.numeric(ran)) stop(paste("'",nm1,"' must be numeric",sep=""))
  if(length(ran)!=2) stop(paste("'",nm1,"' must have length 2",sep=""))
  if(ran[2]<=ran[1]) stop(paste("'",nm1,"[1]' must be < '",nm1,"[2]'",sep=""))
  if(any(vals<ran[1])||any(vals>ran[2])) stop(paste("all time values must lie within interval '",nm1,"'",sep=""))
  return(ran)
}

checktsel <- function(tsel){
  if(!is.vector(tsel)) stop(paste("'tselect' must be a vector",sep=""))
  if(!is.numeric(tsel)) stop(paste("'tselect' must be numeric",sep=""))
  if(length(tsel)==2){
    if(tsel[2]<tsel[1]) stop(paste("'tselect[1]' must be <= 'tselect[2]'",sep=""))
  } else if(length(tsel)==1){
    tsel <- rep(tsel,2)
  } else {
    stop(paste("'tselect' must be of length 1 or 2"))
  }
  return(tsel)
}
  
checktt <- function(tt){
  if(is.null(tt)) stop("observation times must be supplied as 'tt' or as a numeric vector of marks to 'pp'")
  if(!is.vector(tt)) stop("'tt' must be a vector")
  if(!is.numeric(tt)) stop("'tt' must be numeric")
  return(tt)
}

checkxy <- function(xy){
  if(!is.list(xy)) stop("'xy' must be a list")
  if(is.null(xy$x)||is.null(xy$y)) stop("'xy' have vector members 'x' and 'y'")
  if(length(xy$x)!=length(xy$y)) stop("'xy' must be a list with numeric vectors 'x' and 'y' of equal length")
  if((!is.numeric(xy$x))||(!is.numeric(xy$y))) stop("'xy' must have numeric vector members 'x' and 'y'")
  return(xy)
}

checkwei <- function(wei,n){
  if(!is.vector(wei)) stop("'weights' must be a vector")
  if(length(wei)!=n) stop("length of 'weights' must match number of observations")
  if(!is.numeric(wei)) stop("'weights' must be numeric")
  if(any(wei<0)) stop("all 'weights' must be nonnegative")
  return(wei)
}