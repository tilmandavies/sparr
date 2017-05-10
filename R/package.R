#' The sparr Package: SPAtial Relative Risk
#' 
#' Provides functions to estimate fixed and adaptive kernel-smoothed relative
#' risk surfaces via the density-ratio method and perform subsequent inference.
#' 
#' @template version
#' @details
#'
#' Kernel smoothing, and the
#' flexibility afforded by this methodology, provides an attractive approach to
#' estimating complex probability density functions. This is particularly of
#' interest when exploring problems in geographical epidemiology, the study of
#' disease dispersion throughout some spatial region, given a population. The
#' so-called `relative risk surface', constructed as a ratio of estimated case
#' to control densities (Bithell, 1990; 1991), describes the variation in the
#' `risk' of the disease, given the underlying at-risk population. This is a
#' technique that has been applied successfully for mainly exploratory purposes
#' in a number of different examples (see for example Sabel et al., 2000;
#' Prince et al., 2001; Wheeler, 2007).
#' 
#' This package provides functions for bivariate kernel density estimation
#' (KDE), implementing both fixed and `variable' or `adaptive' (Abramson, 1982)
#' smoothing parameter options (see the function documentation for more
#' information). A selection of bandwidth calculators for bivariate KDE and the
#' relative risk function are provided, including one based on the maximal
#' smoothing principle (Terrell, 1990), and others involving a leave-one-out
#' least-squares cross-validation (see below). In addition, the ability to
#' construct asymptotically derived \emph{p}-value surfaces (`tolerance'
#' contours of which signal statistically significant sub-regions of extremity
#' in a risk surface - Hazelton and Davies, 2009; Davies and Hazelton, 2010) as
#' well as some visualisation tools are provided.
#' 
#' The content of \code{sparr} can be broken up as follows:\cr
#' 
#' \emph{Datasets}
#' 
#' \code{\link{pbc}} a case/control planar point pattern
#' (\code{\link[spatstat]{ppp.object}}) concerning liver disease in northern
#' England. Also available are a number of built-in datasets in the
#' \code{\link[spatstat]{spatstat}} package such as
#' \code{\link[spatstat]{chorley}}, which concerns the distribution of
#' laryngeal cancer in an area of Lancashire, England.\cr
#' 
#' \emph{Bandwidth calculators}
#' 
#' \code{\link{OS}} estimation of an isotropic
#' smoothing parameter for fixed-bandwidth bivariate KDE, based on the
#' oversmoothing principle introduced by Terrell (1990).
#' 
#' \code{\link{NS}}
#' estimation of an isotropic smoothing parameter for fixed-bandwidth bivariate
#' KDE, based on the asymptotically optimal value for a normal density
#' (bivariate normal scale rule - see e.g. Wand and Jones, 1995).
#' 
#' \code{\link{LSCV.density}} a least-squares cross-validated (LSCV) estimate
#' of an isotropic fixed bandwidth for bivariate KDE (see e.g. Bowman and
#' Azzalini, 1997).
#' 
#' \code{\link{LSCV.risk}} Estimation of a jointly optimal,
#' common isotropic case-control fixed bandwidth for the kernel-smoothed risk
#' function based on the mean integrated squared error (MISE), a weighted MISE,
#' or the asymptotic MISE (see Kelsall and Diggle, 1995a; Hazelton, 2008;
#' Davies, 2013).
#' 
#' --More in development--
#' 
#' \emph{Spatial functions}
#' 
#' \code{\link{bivariate.density}} kernel density
#' estimate of bivariate data; fixed or adaptive smoothing
#' 
#' \code{\link{multiscale.density}} multi-scale adaptive kernel density
#' estimates for multiple global bandwidths as per Davies and Baddeley
#' (2017).
#' 
#' \code{\link{multiscale.slice}} a single adaptive kernel estimate
#' based on taking a slice from a multi-scale estimate.\cr
#' 
#' \emph{Spatiotemporal functions}
#' 
#' --In development--
#' 
#' \emph{Relative risk and \emph{p}-value surfaces}
#' 
#' \code{\link{risk}}
#' estimation of a (log) spatial relative risk function, either from data or
#' pre-existing bivariate density estimates
#' 
#' \code{\link{tolerance}}
#' calculation of asymptotic or Monte-Carlo \emph{p}-value surfaces
#' 
#' \emph{Printing and summarising objects}
#' 
#' \code{S3} methods (\code{\link{print.bivden}}, \code{\link{print.rrs}},
#' \code{\link{print.msden}}, \code{\link{summary.bivden}},
#' \code{\link{summary.rrs}}, and \code{\link{print.msden}}) are available for
#' the bivariate density, risk function, and multi-scale density objects.
#' 
#' \emph{Visualisation}
#' 
#' \code{S3} methods of the \code{plot} function; see
#' \code{\link{plot.bivden}} for visualising a single bivariate density
#' estimate from \code{\link{bivariate.density}}, \code{\link{plot.rrs}} for
#' visualisation of an estimated relative risk function from
#' \code{\link{risk}}, and \code{\link{plot.msden}} for viewing animations of
#' multi-scale density estimates from \code{\link{multiscale.density}}.
#' 
#' \code{\link{tol.contour}} provides more flexibility for plotting and
#' superimposing tolerance contours upon an existing plot, given output from
#' \code{\link{tolerance}}
#' 
#' @name sparr-package
#' @aliases sparr-package sparr
#' @docType package
#' @section Dependencies: The \code{sparr} package depends upon
#' \code{\link[spatstat]{spatstat}}. In particular, the user should familiarise
#' themselves with \code{\link[spatstat]{ppp}} objects and
#' \code{\link[spatstat]{im}} objects, which are used throughout. For the
#' experimental capabilities involving parallel processing, \code{sparr} also
#' currently imports \code{\link[doParallel]{doParallel}},
#' \code{\link[parallel]{parallel}}, and \code{\link[foreach]{foreach}}.
#' 
#' @author T.M. Davies\cr Dept. of Mathematics & Statistics, University of
#' Otago, Dunedin, New Zealand.\cr
#' J.C. Marshall\cr
#' Institute of Fundamantal Sciences, Massey University, Palmerston North, New Zealand.\cr
#' 
#' Maintainer: T.M.D. \email{tdavies@@maths.otago.ac.nz}
#' 
#' @references
#' Abramson, I. (1982), On bandwidth variation in kernel estimates
#' --- a square root law, \emph{Annals of Statistics}, \bold{10}(4),
#' 1217-1223.
#' 
#' Adler, D. and Murdoch, D. (2009), rgl: 3D visualization device
#' system (OpenGL). R package version 0.87; URL:
#' http://CRAN.R-project.org/package=rgl
#' 
#' Baddeley, A. and Turner, R. (2005),
#' Spatstat: an R package for analyzing spatial point patterns, \emph{Journal
#' of Statistical Software}, \bold{12}(6), 1-42.
#' 
#' Bithell, J.F. (1990), An
#' application of density estimation to geographical epidemiology,
#' \emph{Statistics in Medicine}, \bold{9}, 691-701.
#' 
#' Bithell, J.F. (1991),
#' Estimation of relative risk function,. \emph{Statistics in Medicine},
#' \bold{10}, 1745-1751.
#' 
#' Bowman, A.W. and Azzalini, A. (1997), \emph{Applied
#' Smoothing Techniques for Data Analysis: The Kernel Approach with S-Plus
#' Illustrations.} Oxford University Press Inc., New York. ISBN
#' 0-19-852396-3.
#' 
#' Davies, T. M. (2013), Jointly optimal bandwidth selection
#' for the planar kernel-smoothed density-ratio, \emph{Spatial and
#' Spatio-temporal Epidemiology}, \bold{5}, 51-65.
#' 
#' Davies, T.M. and Baddeley
#' A. (2017), Fast computation of spatially adaptive kernel estimates,
#' \emph{Submitted}.
#' 
#' Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel
#' estimation of spatial relative risk, \emph{Statistics in Medicine},
#' \bold{29}(23) 2423-2437.
#' 
#' Davies, T.M., Hazelton, M.L. and Marshall, J.C.
#' (2011), \code{sparr}: Analyzing spatial relative risk using fixed and
#' adaptive kernel density estimation in \code{R}, \emph{Journal of Statistical
#' Software} \bold{39}(1), 1-14.
#' 
#' Davies, T.M., Jones, K. and Hazelton, M.L.
#' (2016), Symmetric adaptive smoothing regimens for estimation of the spatial
#' relative risk function, \emph{Computational Statistics & Data Analysis},
#' \bold{101}, 12-28.
#' 
#' Hazelton, M. L. (2008), Letter to the editor: Kernel
#' estimation of risk surfaces without the need for edge correction,
#' \emph{Statistics in Medicine}, \bold{27}, 2269-2272.
#' 
#' Hazelton, M.L. and
#' Davies, T.M. (2009), Inference based on kernel estimates of the relative
#' risk function in geographical epidemiology, \emph{Biometrical Journal},
#' \bold{51}(1), 98-109.
#' 
#' Kelsall, J.E. and Diggle, P.J. (1995a), Kernel
#' estimation of relative risk, \emph{Bernoulli}, \bold{1}, 3-16.
#' 
#' Kelsall, J.E. and Diggle, P.J. (1995b), Non-parametric estimation of spatial
#' variation in relative risk, \emph{Statistics in Medicine}, \bold{14},
#' 2335-2342.
#' 
#' Prince, M. I., Chetwynd, A., Diggle, P. J., Jarner, M.,
#' Metcalf, J. V. and James, O. F. W. (2001), The geographical distribution of
#' primary biliary cirrhosis in a well-defined cohort, \emph{Hepatology}
#' \bold{34}, 1083-1088.
#' 
#' Sabel, C. E., Gatrell, A. C., Loytonenc, M.,
#' Maasiltad, P. and Jokelainene, M. (2000), Modelling exposure opportunitites:
#' estimating relative risk for motor disease in Finland, \emph{Social Science
#' & Medicine} \bold{50}, 1121-1137.
#' 
#' Terrell, G.R. (1990), The maximal
#' smoothing principle in density estimation, \emph{Journal of the American
#' Statistical Association}, \bold{85}, 470-477.
#' 
#' Venables, W. N. and Ripley,
#' B. D. (2002). \emph{Modern Applied Statistics with S}, Fourth Edition,
#' Springer, New York.
#' 
#' Wand, M.P. and Jones, C.M., 1995. \emph{Kernel
#' Smoothing}, Chapman & Hall, London.
#' 
#' Wheeler, D. C. (2007), A comparison
#' of spatial clustering and cluster detection techniques for childhood
#' leukemia incidence in Ohio, 1996-2003, \emph{International Journal of Health
#' Geographics}, \bold{6}(13).
#' 
#' @keywords package
NULL

#' @importFrom utils setTxtProgressBar txtProgressBar packageVersion
#' @importFrom stats IQR density dnorm fft optimise pnorm quantile rnorm sd var
#' @importFrom graphics axis box contour pairs par plot points title
#' @importFrom grDevices dev.hold dev.flush
#' @importFrom spatstat.utils prange tapplysum inside.range
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @importFrom ks kde
#' @import spatstat
NULL