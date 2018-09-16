#' Bivariate kernel density/intensity estimation
#' 
#' Provides an isotropic adaptive or fixed bandwidth kernel density/intensity
#' estimate of bivariate/planar/2D data.
#' 
#' Given a data set \eqn{x_1,\dots,x_n} in 2D, the isotropic kernel estimate of
#' its probability density function, \eqn{\hat{f}(x)}{\hat{f}(x)}, is given by
#' \deqn{\hat{f}(y)=n^{-1}\sum_{i=1}^{n}h(x_i)^{-2}K((y-x_i)/h(x_i)) }
#' where \eqn{h(x)}{h(x)} is the bandwidth function, and \eqn{K(.)} is the
#' bivariate standard normal smoothing kernel. Edge-correction factors (not
#' shown above) are also implemented.
#' 
#' \describe{
#'   \item{\bold{Fixed}}{
#'     The classic fixed bandwidth kernel estimator is used when
#'     \code{adapt = FALSE}. This amounts to setting \eqn{h(u)=}\code{h0} for all \eqn{u}.
#'     Further details can be found in the documentation for \code{\link[spatstat]{density.ppp}}.}
#'   \item{\bold{Adaptive}}{Setting \code{adapt = TRUE} requests computation of Abramson's (1982)
#'     variable-bandwidth estimator. Under this framework, we have
#'     \eqn{h(u)=}\code{h0}*min[\eqn{\tilde{f}(u)^{-1/2}},\eqn{G*}\code{trim}]/\eqn{\gamma},
#'     where \eqn{\tilde{f}(u)} is a fixed-bandwidth kernel density estimate
#'     computed using the pilot bandwidth \code{hp}.
#'     \itemize{
#'       \item Global smoothing of the variable bandwidths is controlled with the global bandwidth
#'         \code{h0}.
#'       \item In the above statement, \eqn{G} is the geometric mean of the
#'         ``bandwidth factors'' \eqn{\tilde{f}(x_i)^{-1/2}}; \eqn{i=1,\dots,n}. By
#'         default, the variable bandwidths are rescaled by \eqn{\gamma=G}, which is
#'         set with \code{gamma.scale = "geometric"}. This allows \code{h0} to be
#'         considered on the same scale as the smoothing parameter in a fixed-bandwidth
#'         estimate i.e. on the scale of the recorded data. You can use any other
#'         rescaling of \code{h0} by setting \code{gamma.scale} to be any scalar
#'         positive numeric value; though note this only affects \eqn{\gamma} -- see
#'         the next bullet. When using a scale-invariant \code{h0}, set
#'         \code{gamma.scale = 1}.
#'      \item The variable bandwidths must be trimmed to
#'         prevent excessive values (Hall and Marron, 1988). This is achieved through
#'         \code{trim}, as can be seen in the equation for \eqn{h(u)} above. The
#'         trimming of the variable bandwidths is universally enforced by the geometric
#'         mean of the bandwidth factors \eqn{G} independent of the choice of
#'         \eqn{\gamma}. By default, the function truncates bandwidth factors at five
#'         times their geometric mean. For stricter trimming, reduce \code{trim}, for
#'         no trimming, set \code{trim = Inf}.
#'      \item For even moderately sized data sets
#'         and evaluation grid \code{resolution}, adaptive kernel estimation can be
#'         rather computationally expensive. The argument \code{davies.baddeley} is
#'         used to approximate an adaptive kernel estimate by a sum of fixed bandwidth
#'         estimates operating on appropriate subsets of \code{pp}. These subsets are
#'         defined by ``bandwidth bins'', which themselves are delineated by a quantile
#'         step value \eqn{0<\delta<1}. E.g. setting \eqn{\delta=0.05} will create 20
#'         bandwidth bins based on the 0.05th quantiles of the Abramson variable
#'         bandwidths. Adaptive edge-correction also utilises the partitioning, with
#'         pixel-wise bandwidth bins defined using the value \eqn{0<\beta<1}, and the
#'         option to decrease the resolution of the edge-correction surface for
#'         computation to a [\eqn{L} \eqn{\times}{x} \eqn{L}] grid, where \eqn{0 <L
#'         \le} \code{resolution}. If \code{davies.baddeley} is supplied as a vector of
#'         length 3, then the values \code{[1], [2], and [3]} correspond to the
#'         parameters \eqn{\delta}, \eqn{\beta}, and \eqn{L_M=L_N} in Davies and
#'         Baddeley (2018). If the argument is simply a single numeric value, it is
#'         used for both \eqn{\delta} and \eqn{\beta}, with
#'         \eqn{L_M=L_N=}\code{resolution} (i.e. no edge-correction surface
#'         coarsening).
#'      \item Computation of leave-one-out values (when
#'         \code{leaveoneout = TRUE}) is done by brute force, and is therefore very
#'         computationally expensive for adaptive smoothing. This is because the
#'         leave-one-out mechanism is applied to both the pilot estimation and the
#'         final estimation stages. Experimental code to do this via parallel
#'         processing using the \code{\link{foreach}} routine is implemented.
#'         Fixed-bandwidth leave-one-out can be performed directly in
#'         \code{\link[spatstat]{density.ppp}}.
#'     }
#'   }
#' }
#' 
#' @aliases bivariate.density bivden
#' 
#' @param pp An object of class \code{\link[spatstat]{ppp}} giving the observed
#'   2D data set to be smoothed.
#' @param h0 Global bandwidth for adaptive smoothing or fixed bandwidth for
#'   constant smoothing. A numeric value > 0.
#' @param hp Pilot bandwidth (scalar, numeric > 0) to be used for fixed
#'   bandwidth estimation of a pilot density in the case of adaptive smoothing.
#'   If \code{NULL} (default), it will take on the value of \code{h0}. Ignored
#'   when \code{adapt = FALSE} or if \code{pilot.density} is supplied as a
#'   pre-defined pixel image.
#' @param adapt Logical value indicating whether to perform adaptive kernel
#'   estimation. See `Details'.
#' @param resolution Numeric value > 0. Resolution of evaluation grid; the
#'   density/intensity will be returned on a [\code{resolution} \eqn{\times}{x}
#'   \code{resolution}] grid.
#' @param gamma.scale Scalar, numeric value > 0; controls rescaling of the
#'   variable bandwidths. Defaults to the geometric mean of the bandwidth factors
#'   given the pilot density (as per Silverman, 1986). See `Details'.
#' @param edge Character string giving the type of edge correction to employ.
#'   \code{"uniform"} (default) corrects based on evaluation grid coordinate and
#'   \code{"diggle"} reweights each observation-specific kernel. Setting
#'   \code{edge = "none"} requests no edge correction. Further details can be
#'   found in the documentation for \code{\link[spatstat]{density.ppp}}.
#' @param weights Optional numeric vector of nonnegative weights corresponding to
#'   each observation in \code{pp}. Must have length equal to \code{npoints(pp)}.
#' @param intensity Logical value indicating whether to return an intensity
#'   estimate (integrates to the sample size over the study region), or a density
#'   estimate (default, integrates to 1).
#' @param trim Numeric value > 0; controls bandwidth truncation for adaptive
#'   estimation. See `Details'.
#' @param xy Optional alternative specification of the evaluation grid; matches
#'   the argument of the same tag in \code{\link[spatstat]{as.mask}}. If
#'   supplied, \code{resolution} is ignored.
#' @param pilot.density An optional pixel image (class
#'   \code{\link[spatstat]{im}}) giving the pilot density to be used for
#'   calculation of the variable bandwidths in adaptive estimation, \bold{or} a
#'   \code{\link[spatstat]{ppp.object}} giving the data upon which to base a
#'   fixed-bandwidth pilot estimate using \code{hp}. If used, the pixel image
#'   \emph{must} be defined over the same domain as the data given
#' \code{resolution} or the supplied pre-set \code{xy} evaluation grid;
#'   \bold{or} the planar point pattern data must be defined with respect to the
#'   same polygonal study region as in \code{pp}.
#' @param leaveoneout Logical value indicating whether to compute and return
#'   the value of the density/intensity at each data point for an adaptive
#'   estimate. See `Details'.
#' @param parallelise Numeric argument to invoke parallel processing, giving
#'   the number of CPU cores to use when \code{leaveoneout = TRUE}. Experimental.
#'   Test your system first using \code{parallel::detectCores()} to identify the
#'   number of cores available to you.
#' @param davies.baddeley An optional numeric vector of length 3 to control
#'   bandwidth partitioning for approximate adaptive estimation, giving the
#'   quantile step values for the variable bandwidths for density/intensity and
#'   edge correction surfaces and the resolution of the edge correction surface.
#'   May also be provided as a single numeric value. See `Details'.
#' @param verbose Logical value indicating whether to print a function progress
#'   bar to the console when \code{adapt = TRUE}.
#' 
#' @return If \code{leaveoneout = FALSE}, an object of class \code{"bivden"}.
#' This is effectively a list with the following components:
#' \item{z}{The
#' resulting density/intensity estimate, a pixel image object of class
#' \code{\link[spatstat]{im}}.}
#' 
#' \item{h0}{A copy of the value of \code{h0}
#' used.} \item{hp}{A copy of the value of \code{hp} used.}
#' 
#' \item{h}{A numeric
#' vector of length equal to the number of data points, giving the bandwidth
#' used for the corresponding observation in \code{pp}.}
#' 
#' \item{him}{A pixel
#' image (class \code{\link[spatstat]{im}}), giving the `hypothetical' Abramson
#' bandwidth at each pixel coordinate conditional upon the observed data.
#' \code{NULL} for fixed-bandwidth estimates.}
#' 
#' \item{q}{Edge-correction
#' weights; a pixel \code{\link[spatstat]{im}}age if \code{edge = "uniform"}, a
#' numeric vector if \code{edge = "diggle"}, and \code{NULL} if \code{edge =
#' "none"}.}
#' 
#' \item{gamma}{The value of \eqn{\gamma} used in scaling the
#' bandwidths. \code{NA} if a fixed bandwidth estimate is computed.}
#' 
#' \item{geometric}{The geometric mean \eqn{G} of the untrimmed bandwidth
#' factors \eqn{\tilde{f}(x_i)^{-1/2}}. \code{NA} if a fixed bandwidth estimate
#' is computed.}
#' 
#' \item{pp}{A copy of the \code{\link[spatstat]{ppp.object}}
#' initially passed to the \code{pp} argument, containing the data that were
#' smoothed.}
#' 
#' 
#' Else, if \code{leaveoneout = TRUE}, simply a numeric vector of length equal to the
#' number of data points, giving the leave-one-out value of the function at the
#' corresponding coordinate.
#' 
#' @author T.M. Davies and J.C. Marshall
#' 
#' @references
#' Abramson, I. (1982). On bandwidth variation in kernel estimates
#' --- a square root law, \emph{Annals of Statistics}, \bold{10}(4),
#' 1217-1223.
#' 
#' Davies, T.M. and Baddeley A. (2018), Fast computation of
#' spatially adaptive kernel estimates, \emph{Statistics and Computing}, \bold{28}(4), 937-956.
#' 
#' Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative
#' risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.
#' 
#' Davies, T.M., Jones, K. and Hazelton, M.L. (2016), Symmetric adaptive smoothing
#' regimens for estimation of the spatial relative risk function,
#' \emph{Computational Statistics & Data Analysis}, \bold{101}, 12-28.
#' 
#' Diggle, P.J. (1985), A kernel method for smoothing point process data,
#' \emph{Journal of the Royal Statistical Society, Series C}, \bold{34}(2),
#' 138-147.
#' 
#' Hall P. and Marron J.S. (1988) Variable window width kernel
#' density estimates of probability densities. \emph{Probability Theory and
#' Related Fields}, \bold{80}, 37-49.
#' 
#' Marshall, J.C. and Hazelton, M.L. (2010) Boundary kernels for adaptive density
#' estimators on regions with irregular boundaries, \emph{Journal of Multivariate
#' Analysis}, \bold{101}, 949-963.
#' 
#' Silverman, B.W. (1986), \emph{Density Estimation for
#' Statistics and Data Analysis}, Chapman & Hall, New York.
#' 
#' Wand, M.P. and Jones, C.M., 1995. \emph{Kernel Smoothing}, Chapman & Hall, London.
#' 
#' @examples
#' 
#' data(chorley) # Chorley-Ribble data from package 'spatstat'
#' 
#' # Fixed bandwidth kernel density; uniform edge correction
#' chden1 <- bivariate.density(chorley,h0=1.5) 
#' 
#' # Fixed bandwidth kernel density; diggle edge correction; coarser resolution
#' chden2 <- bivariate.density(chorley,h0=1.5,edge="diggle",resolution=64) 
#' 
#' \donttest{
#' # Adaptive smoothing; uniform edge correction
#' chden3 <- bivariate.density(chorley,h0=1.5,hp=1,adapt=TRUE)
#' 
#' # Adaptive smoothing; uniform edge correction; partitioning approximation
#' chden4 <- bivariate.density(chorley,h0=1.5,hp=1,adapt=TRUE,davies.baddeley=0.025)
#'  
#' par(mfrow=c(2,2))
#' plot(chden1);plot(chden2);plot(chden3);plot(chden4)  
#' }
#' 
#' @export
bivariate.density <- function(pp,h0,hp=NULL,adapt=FALSE,resolution=128,gamma.scale="geometric",edge=c("uniform","diggle","none"),weights=NULL,intensity=FALSE,trim=5,xy=NULL,pilot.density=NULL,leaveoneout=FALSE,parallelise=NULL,davies.baddeley=NULL,verbose=TRUE){
	if(!inherits(pp,"ppp")) stop("data argument 'pp' must be of spatstat class \"ppp\"; see ?ppp")

	W <- Window(pp)
	if(!is.null(xy)){
	  xy <- checkxy(xy)
	  dimyx <- NULL
	  resolution <- length(xy$x)
	} else {
	  resolution <- checkit(resolution,"'resolution'")
	  dimyx <- rep(resolution,2)
	}

	edge <- checkedge(edge)
	h0 <- checkit(h0,"'h0'")
	if(!is.null(trim)&&!is.na(trim)) trim <- checkit(trim,"'trim'")
	
	n <- npoints(pp)
  if(!is.null(weights)) weights <- checkwei(weights,n)
	
	if(adapt){
		if(is.null(hp)) hp <- h0
		else hp <- checkit(hp,"'hp'")
		
		if(leaveoneout) return(bivden.LOO(pp,h0,hp,(edge=="uniform"||edge=="diggle"),gamma.scale,trim,resolution,parallelise,weights,0)[[1]])
    
		pd <- pilot.density
		pilot.data <- pp
		if(!is.null(pd)){
			if(is.im(pd)){
				if(is.null(xy)){
			    		if(!all(dim(pd)==resolution)) stop("'pilot.density' image resolution must strictly have 'resolution' x 'resolution' pixels")
				} else {
			    		if((!all(pd$xcol==xy$x))||(!all(pd$yrow==xy$y))) stop("'pilot.density' xcol and yrow must strictly match coords in 'xy'")
				}
		  	pilot.density[pd<=0] <- min(pd[pd>0])
		  	hp <- NULL
			} else if(is.ppp(pd)){
				pilot.data <- pd
				if(!identical_windows(Window(pp),Window(pilot.data))) stop("'pilot.density' window must be identical to 'pp' window")
				pilot.density <- density(pilot.data,sigma=hp,edge=(edge=="uniform"||edge=="diggle"),diggle=(edge=="diggle"),dimyx=dimyx,xy=xy,positive=TRUE)
			} else {
				stop("'pilot.density' must be an object of class \"im\" or \"ppp\"")
			}
		} else {
			pilot.density <- density(pp,sigma=hp,edge=(edge=="uniform"||edge=="diggle"),diggle=(edge=="diggle"),dimyx=dimyx,xy=xy,positive=TRUE,weights=weights)
		}
		
		pilot.density.spec <- safelookup(pilot.density,pp,warn=FALSE)
		pi.int <- integral(pilot.density)
		pilot.density <- pilot.density/pi.int
		pilot.density.spec <- pilot.density.spec/pi.int
		pspec <- pilot.density.spec^(-0.5)
		gamma <- processgamma(gamma.scale,safelookup(pilot.density,pilot.data,warn=FALSE)) #'was: processgamma(gamma.scale,pilot.density.spec)
		gs <- gspd <- exp(mean(log(pspec))) 
		if(!is.null(pd)) gspd <- exp(mean(log(safelookup(pilot.density,pilot.data,warn=FALSE)^(-0.5))))
		
		# PREVIOUS TRIMMING REGIMEN #
		# h.spec <- h0*pilot.density.spec^(-0.5)/gamma
		# h.hypo <- h0*pilot.density^(-0.5)/gamma
		# if(is.null(trim)) beta.h <- 5*median(h.spec,na.rm=TRUE)
		# else if(is.na(trim)) beta.h <- max(h.hypo,na.rm=TRUE)
		# else beta.h <- trim
		# h.spec[h.spec>beta.h] <- beta.h
		# h.hypo[h.hypo>beta.h] <- beta.h
		# h.hypo.mat <- as.matrix(h.hypo)
		
		# NEW TRIMMING REGIMEN #
		# h.spec <- h0*pmin(pspec/gamma,trim)  ### Generalised below for numeric gamma argument vals. Trimming is universally determined by the geometric mean 'gs', regardless of 'gamma.scale' ###
		# h.hypo <- h0*im(matrix(pmin(as.vector(as.matrix(pilot.density^(-0.5)))/gamma,trim),resolution,resolution),xcol=pilot.density$xcol,yrow=pilot.density$yrow)
		
		h.spec <- h0*pmin(pspec,trim*gspd)/gamma
		h.hypo <- h0*im(matrix(pmin(as.vector(as.matrix(pilot.density^(-0.5))),trim*gspd),resolution,resolution)/gamma,xcol=pilot.density$xcol,yrow=pilot.density$yrow)
		h.hypo.mat <- as.matrix(h.hypo)
		h.hypo.vec <- as.numeric(t(h.hypo.mat))

		if(!is.null(davies.baddeley)){
			db <- checkdb(davies.baddeley)
			db.result <- adens(x=pp,bwim=h.hypo,bwpts=h.spec,resolution=resolution,intensity=intensity,edge=(edge=="uniform"||edge=="diggle"),diggle=(edge=="diggle"),weights=weights,hstep=db[1],qstep=db[2],qres=ifelse(is.na(db[3]),resolution,db[3]),verbose)
			result <- list(z=db.result$result,h0=h0,hp=hp,h=h.spec,him=h.hypo,q=db.result$edg,gamma=gamma,geometric=gs,pp=pp)
			class(result) <- "bivden"		
			return(result)
		}

		evalxy <- as.matrix(expand.grid(pilot.density$xcol,pilot.density$yrow))
		notin <- !inside.owin(x=evalxy[,1],y=evalxy[,2],w=W)
		surf <- rep(NA,nrow(evalxy))
		ef <- NULL

		# we need only evaluate on points inside the owin
		evalxy.in <- evalxy[!notin,]
		h.hypo.in <- h.hypo.vec[!notin]
		surf.in   <- numeric(nrow(evalxy.in))
    
		if(is.null(weights)) weights <- rep(1,n)
		if(edge=="uniform"){
			qhz.in <- numeric(nrow(evalxy.in))
			if(verbose) pb <- txtProgressBar(0,nrow(evalxy.in))
			for(i in 1:nrow(evalxy.in)){
				gxy <- kernel2d(evalxy.in[,1]-evalxy.in[i,1], evalxy.in[,2]-evalxy.in[i,2], h.hypo.in[i])
				qhz.in[i] <- dintegral(gxy, pilot.density$xstep, pilot.density$ystep)
		    ivals <- kernel2d(pp$x-evalxy.in[i,1], pp$y-evalxy.in[i,2], h.spec)

        # if(!intensity) surf.in[i] <- mean(ivals)/qhz.in[i]
        # else 
        surf.in[i] <- sum(weights*ivals)/qhz.in[i]
        if(verbose) setTxtProgressBar(pb,i)
			}
			if(verbose) close(pb)
			qhz <- rep(NA,resolution^2)
			qhz[!notin] <- qhz.in
			ef <- im(matrix(qhz,resolution,resolution,byrow=TRUE),xcol=pilot.density$xcol,yrow=pilot.density$yrow)
		}
		
		if(edge=="diggle"){
		  qx <- rep(1,n)
		  if(verbose) pb <- txtProgressBar(0,n+nrow(evalxy.in))
		  for(i in 1:n){
		    pxy <- kernel2d(evalxy.in[,1]-pp$x[i], evalxy.in[,2]-pp$y[i], h.spec[i])
		    qx[i] <- dintegral(pxy, pilot.density$xstep, pilot.density$ystep)
		    if(verbose) setTxtProgressBar(pb,i)
		  }
		  
		  for(i in 1:nrow(evalxy.in)){
		    ivals <- kernel2d(pp$x-evalxy.in[i,1], pp$y-evalxy.in[i,2], h.spec)
		    # if(!intensity) surf.in[i] <- mean(ivals/qx)
		    # else 
		    surf.in[i] <- sum(weights*ivals/qx)
		    if(verbose) setTxtProgressBar(pb,n+i)
		  }
		  if(verbose) close(pb)
		  ef <- qx
		}
		  
		if(edge=="none"){
		  if(verbose) pb <- txtProgressBar(0,nrow(evalxy.in))
		  for(i in 1:nrow(evalxy.in)){
		    ivals <- kernel2d(pp$x-evalxy.in[i,1], pp$y-evalxy.in[i,2], h.spec)
		    # if(!intensity) surf.in[i] <- mean(ivals)
		    # else 
		    surf.in[i] <- sum(weights*ivals)
		    if(verbose) setTxtProgressBar(pb,i)
		  }
		  if(verbose) close(pb)
		}
		
		surf[!notin] <- surf.in
		surf <- im(matrix(surf,resolution,resolution,byrow=TRUE),xcol=pilot.density$xcol,yrow=pilot.density$yrow)

	} else {
    h.spec <- rep(h0,n)
	  h.hypo <- NULL
	  gs <- gamma <- NA
    
	  dens <- density.ppp(pp,sigma=h0,dimyx=dimyx,xy=xy,edge=(edge=="diggle"||edge=="uniform"),diggle=(edge=="diggle"),weights=weights,spill=1)
		surf <- dens$raw[W,drop=FALSE]
		ef <- dens$edg[W,drop=FALSE]
		ef[ef>1] <- 1
		
		if(edge=="diggle"){
		  ef <- safelookup(ef,pp,warn=FALSE)
		} else if(edge=="uniform"){
		  surf <- surf/ef
		  surf[surf<0] <- 0
		} else {
		  ef <- NULL
		}
	}
	
	if(!intensity) surf <- surf/integral(surf)
	result <- list(z=surf,h0=h0,hp=hp,h=h.spec,him=h.hypo,q=ef,gamma=gamma,geometric=gs,pp=pp)
	class(result) <- "bivden"		
	
	return(result)
}
