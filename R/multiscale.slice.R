multiscale.slice <- function(msob,h0,checkargs=TRUE){
  if(checkargs){
    if(!inherits(msob,"msden")) stop("'msob' must be of class \"msden\"")
    h0 <- checkit(h0,"'h0'")
    aran <- msob$h0range
    if(!inside.range(h0,aran)) stop(paste("requested 'h0' outside available range of",prange(aran)))
  }
  
  available.h0 <- as.numeric(names(msob$z))
  zz <- msob$z
  hh <- msob$him
  qq <- msob$q
  
  if(any(available.h0==h0)){
    index <- which(available.h0==h0)
    zres <- zz[[index]]
    hres <- hh[[index]]
    qres <- qq[[index]]
  } else {
    marker <- which(available.h0>h0)[1]
    mindex <- c(marker-1,marker)
    hint <- available.h0[mindex]
    move <- (h0-hint[1])/diff(hint)
    zdiff <- zz[[mindex[2]]]-zz[[mindex[1]]]
    hdiff <- hh[[mindex[2]]]-hh[[mindex[1]]]
    qdiff <- qq[[mindex[2]]]-qq[[mindex[1]]]
    zres <- zz[[mindex[1]]]+move*zdiff
    hres <- hh[[mindex[1]]]+move*hdiff
    if(!is.null(qq)) qres <- qq[[mindex[1]]]+move*qdiff
    else qres <- NULL
  }
  
  result <- list(z=zres,h0=h0,hp=msob$hp,h=msob$h/msob$h0*h0,him=hres,q=qres,gamma=msob$gamma,geometric=msob$geometric,pp=msob$pp,fromms=TRUE)
  class(result) <- "bivden"
  return(result)
}