processgamma <- function(gamma.scale,pds){
  if(length(gamma.scale)==1){
    if(is.numeric(gamma.scale)&&gamma.scale>0){
      gamma <- gamma.scale
    } else if(gamma.scale=="geometric"){
      gamma <- exp(mean(log(pds^(-0.5))))
    } else {
      gamma <- 1
      warning("invalid 'gamma.scale' value -- assuming gamma=1")
    }
  } else if(is.null(gamma.scale)){
    gamma <- 1
  } else {
    gamma <- 1
    warning("invalid 'gamma.scale' value -- assuming gamma=1")
  }
  return(gamma)
}