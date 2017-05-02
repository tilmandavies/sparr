checkran <- function(ran,nm){
  if(!is.vector(ran)) stop(paste("'",nm,"' must be a vector",sep=""))
  if(!is.numeric(ran)) stop(paste("'",nm,"' must be numeric",sep=""))
  if(length(ran)!=2) stop(paste("'",nm,"' must have length 2",sep=""))
  if(ran[2]<=ran[1]) stop(paste("'",nm,"[1]' must be < '",nm,"[2]'",sep=""))
  if(any(ran<=0)) stop(paste("'",nm,"' must be wholly positive",sep=""))
  return(ran)
}