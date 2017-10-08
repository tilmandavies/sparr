# sparr 2.1-13

* Minor documentation link corrections

# sparr 2.1-12

* Includes `weights` argument for `bivariate.density`

# sparr 2.1-11

* Citation information and NEWS file changes.

# sparr 2.1-10

* MAJOR CHANGES since versions <= 0.3-8; NO BACKWARDS COMPATIBILITY

* Accompanying tutorial currently submitted for publication; 
  contact maintainer for preprint or see https://arxiv.org/abs/1707.06888

* User is directed to the examples as part of `help("sparr")` for further assistance
      
## New Functions

* `LIK.density`, `BOOT.density`: New functions for bivariate density bandwidth selection.

* `multiscale.density`, `multiscale.slice`: Multi-scale adaptive kernel density estimation.

* `OS.spattemp`, `NS.spattemp`, `LSCV.spattemp`, `LIK.spattemp`, `BOOT.spattemp`: Bandwidth selection for spatiotemporal bandwidth selection.

* `spattemp.density`, `spattemp.risk`, `spattemp.slice`: Spatiotemporal density and relative risk estimation.
    
* `plot.stden`, `plot.msden`, `plot.rrst`, `tol.contour`: New plotting functions.
    
* `burk`, `fmd`: New datasets.

## Significant user-visible changes

* `bivariate.density`, `risk`, `tolerance`, `LSCV.density`, `LSCV.risk`, `OS`, `NS`, `plot` usage changed.

* Significant speed improvements.
