tol.contour <- function(pim, test = c("upper", "lower", "two-sided"), ...){
  if(!inherits(pim,"im")){
    stop("'pim' must be an object of class 'im', typically arising from a call to 'tolerance2'")
  }
  tt <- t(as.matrix(pim))
  test <- test[1]
  if(test=="lower"){
    tt <- 1-tt
  } else if(test=="two-sided"){
    tt <- 2*pmin(tt,1-tt)
  } else if(test!="upper"){
    stop("invalid 'test'")
  }
  contour(x=pim$xcol,y=pim$yrow,z=tt,...)
}