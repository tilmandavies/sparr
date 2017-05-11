#' 2D Fast fourier wrapper around fftwtools or stats package
#' 
#' Utilises the Fastest fourier transform in the West via the fftwtools
#' package if available
#' 
#' @param x the data to transform
#' @param inverse whether it should compute the inverse transform (defaults to FALSE)
#' @param fftw whether fftwtools package is available.
#' @return the fast fourier (inverse) transform of the same size as x, but complex.
fft2d <- function(x, inverse=FALSE, fftw = fftw_available()) {
  if (fftw) {
    fftwtools::fftw2d(data=x, inverse=inverse)
  } else {
    stats::fft(z=x, inverse=inverse)
  }
}

fftw_available <- function() {
  requireNamespace("fftwtools", quietly=TRUE)
}
