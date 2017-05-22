# Bivariate kernel function centered at 0,0 with isotropic bandwidth h
kernel2d <- function(x, y, h) {
  1/(2*pi*h*h)*exp(-0.5*(x*x+y*y)/(h*h))
}

# Approximation of definite integral of f over a grid with gridsize dx x dy via summation
dintegral <- function(f, dx, dy) {
  sum(f)*dx*dy
}

# Fast fourier transform of a 2D Gaussian based on the continuous FT, where
# sigma is the standard deviation, ds,dt is the time resolution in x,y
# and 2n*2n the number of points

kernel2d_fft <- function(sigma, ds, dt, n) {
  kernel_fft <- function(sigma., dt., n.) {
    f <- c(0:(n.-1),-n.:-1)*pi*sigma./(n.*dt.)
    exp(-0.5*f*f)/dt.
  }
  fZ.x <- kernel_fft(sigma, ds, n)
  fZ.y <- kernel_fft(sigma, dt, n)
  fZ.y %*% t(fZ.x)
}