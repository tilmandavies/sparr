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