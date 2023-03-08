# sparr 2.3-10

* Incorporated new shrinkage estimator of spatial relative risk

## Significant user-visible changes

* `risk`: New arguments `shrink` and `shrink.args` implemented for shrinkage estimator

# sparr 2.2-17

* Further further updated to handle new `spatstat` structure

# sparr 2.2-16

* Further updated to handle new `spatstat` structure

# sparr 2.2-15

* Updated internals to handle new `spatstat` structure

# sparr 2.2-14

* Improved handling of `rmdiag` in `BOOT.density` for large datasets
* Increased text output when running `BOOT.density` for fixed bandwidth with `auto.optim=TRUE`

# sparr 2.2-13

* Changed default of `hold` argument in `SLIK.adapt`
* Finer control when plotting `rrst` objects
* Reference updates

## Significant user-visible changes

* `SLIK.adapt`: Argument `hold` now defaults to `TRUE`
* `plot.rrst`: New argument `expscale` for raw-risk scale on image plots of log-relative risk surfaces 

# sparr 2.2-12

* Added `optim.control` argument to `SLIK.adapt`
* Added internal maximum bandwidth limit to `SLIK.adapt`
* Minor bugs squashed, doc updates

# sparr 2.2-11

* Updated citation information and package description

# sparr 2.2-01

* New experimental function `SLIK.adapt` for simultaneous global-pilot bandwidth selection for adaptive densities
* Newly visible function `rimpoly` for random point generation
* Additional control over small bandwidth action in leave-one-out bandwidth selectors
* Improved edge-correction calculations for brute force adaptive leave-one-out
* Updated citation information and documentation corrections

## New Functions

* `rimpoly`: Random spatial point generation based on a pixel image, returned with a polygonal window
* `SLIK.adapt`: Simultaneous bandwidth selector for global and pilot bandwidth for adaptive density estimates based on likelihood cross-validation

## Significant user-visible changes

* `LIK.density`, `LSCV.density`: New argument `zero.action` to provide greater control over the behaviour of leave-one-out calculations at very small bandwidths.

# sparr 2.1-14

* Updated citation information

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
