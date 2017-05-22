# Bivariate kernel function centered at 0,0 with isotropic bandwidth h
kernel2d <- function(x, y, h) {
  1/(2*pi*h*h)*exp(-0.5*(x*x+y*y)/(h*h))
}

# Approximation of definite integral of f over a grid with gridsize dx x dy via summation
dintegral <- function(f, dx, dy) {
  sum(f)*dx*dy
}