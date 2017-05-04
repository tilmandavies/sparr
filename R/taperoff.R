
taperoff <- function(x, zeropoint=0, onepoint=1,
                     type=c("smooth", "cosine")) {
  type <- match.arg(type)
  #   cosine taper is standard in engineering (apparently)
  #
  #   smooth taper is the pure mathematicians' favorite example
  #                of a Smooth Partition of Unity
  y <- (x-zeropoint)/(onepoint - zeropoint)
  z <- switch(type,
              cosine = ifelse(y <= 0, 0,
                              ifelse(y >= 1, 1,
                                     (1 - cos(pi * y))/2)),
              smooth = ifelse(y <= 0, 0,
                              ifelse(y >= 1, 1,
                                     exp(-1/y)/(exp(-1/y) + exp(-1/(1-y))))))
  return(z)
}
