#' 2D fast-Fourier wrapper around 'fftwtools' or 'stats' package
#' 
#' Utilises the Fastest Fourier Transform in the West (FFTW) via the 'fftwtools'
#' package if available, else reverts to built-in functionality
#' 
#' This function is called wherever \code{sparr} seeks to perform a 2D fast-Fourier
#' transform. Where available, computational expense is noticeably reduced by appealing to routines 
#' in the independent `FFTW' toolbox. The user is encouraged to install the corresponding R package \code{fftwtools} from CRAN;
#' this function will automatically detect and use the faster option, otherwise will 
#' defer to the built-in \code{\link[stats]{fft}}.
#' 
#' 
#' @param x A numeric matrix to be transformed.
#' @param inverse Whether it should compute the inverse transform (defaults to \code{FALSE}).
#' @param fftw Whether the \code{fftwtools} R package is available.
#' @return The fast-Fourier (inverse) transform. A complex-valued matrix of the same size as \code{x}.
#' 
#' @author J.C. Marshall
#' 
#' @examples
#' \donttest{
#' 
#' # System check
#' sparr:::fftw_available()
#' 
#' system.time(fft(matrix(1:2000^2,2000)))
#' system.time(fft2d(matrix(1:2000^2,2000)))
#' }
#' 
#' @export
fft2d <- function(x, inverse=FALSE, fftw = sparr:::fftw_available()) {
  if (fftw) {
    fftwtools::fftw2d(data=x, inverse=inverse)
  } else {
    stats::fft(z=x, inverse=inverse)
  }
}

fftw_available <- function() {
  yeahnah <- requireNamespace("fftwtools", quietly=TRUE)
  if(yeahnah) yeahnah <- packageVersion("fftwtools") >= "0.9-8"
  return(yeahnah)
}

# The user is encouraged to install the `FFTW' tools 
# from \url{http://www.fftw.org} and the corresponding R package \code{fftwtools} from CRAN;
# this function will automatically detect and use the faster option, otherwise will 
# defers to the built-in \code{\link[stats]{fft}}.
#