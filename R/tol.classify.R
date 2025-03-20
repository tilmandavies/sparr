#' Classification by \emph{p}-value surfaces
#' 
#' Classifies observed case/control points according to an estimated 
#' \emph{p}-value surface.
#' 
#' 
#' This function takes in a relative risk surface computed with
#'  \code{\link{risk}} and corresponding p-value surface (the latter used for
#'  drawing tolerance contours), and attempts to classify both the case and 
#'  control points as either falling within or without contours drawn at a level
#'  of \code{cutoff}. Points that fall 'inside' the contours are deemed to be
#'  associated with p-values less than or equal to \code{cutoff} and hence are
#'  usually interpreted as being in spatial areas of significant risk. This is
#'  useful for identifying characteristics of points that fall inside
#'  'pockets of significance' as delineated by tolerance contours.
#'  
#'  Upon execution, the function first inspects the \code{rs} object to
#'  determine whether it possesses a \code{P} component (i.e. an internally
#'  computed p-value surface provided when \code{\link{risk}} is called with 
#'  optional argument \code{tolerate=TRUE}). If it exists, this is used. If not,
#'  the function then looks to see if the \code{pim} argument has been supplied.
#'  If it has, it must be a pixel \code{\link[spatstat.geom]{im}}age compatible
#'  with the risk surface in \code{rs$rr}. If neither \code{rs$P} nor \code{pim}
#'  is present, the function internally calls \code{\link{tolerance}} with 
#'  arguments supplied to \code{...} to produce the desired surface.
#'  
#'  The return object is a list that splits each of the case and control \code{\link[spatstat.geom]{ppp}}
#'  data objects (these are stored as \code{rs$f$pp} and \code{rs$g$pp}) in the
#'  originally supplied risk surface object) into two constituent \code{\link[spatstat.geom]{ppp}}
#'  objects -- one comprising the points inside the \code{cutoff} contours (\code{fin} and \code{gin}), the 
#'  other for those points outside the \code{cutoff} contours (\code{fout} and \code{gout}).
#'  In addition, the index values of the original data objects \code{rs$f$pp} and 
#'  \code{rs$g$pp} that correspond to the points in \code{fin} and \code{gin} are
#'  provided as numeric vectors (\code{findex} and \code{gindex}). These objects
#'  are useful if you need to cross-reference data-specific characteristics from 
#'  some other (corresponding) data set.
#'  
#'  Further supplied in the returned list are quantities describing the overall classification
#'  structure (\code{pcmask}), as well as contour-specific identification and classification
#'  (\code{finsplit}, \code{ginsplit}, \code{pcpolys}). The \code{pcpolys} object can be plotted
#'  to illustrate the unique contour IDs with \code{\link{tol.classplot}}.
#'  
#' 
#' @param rs An object of class \code{\link{rrs}} giving the estimated relative
#'   risk function of the case-control points to be classified. 
#' @param cutoff A numeric value between 0 and 1, defining the cutoff p-value
#'   used to classify points; defaults to 0.05.
#' @param pim A pixel \code{\link[spatstat.geom]{im}}age defining the p-value 
#'   surface with respect to which the observations are to be classified. 
#'   Typically the result of a call to \code{\link{tolerance}}. Ignored if 
#'   \code{rs} possesses a non-\code{NULL} \code{P} component (which takes
#'   precedence).
#' @param ...  Arguments to be passed to \code{\link{tolerance}} in order to
#'   compute the desired p-value surface in the event neither \code{rs$P} nor
#'   \code{pim} exist. Ignored otherwise.
#'   
#' @return A list of ten components:
#' 
#' \item{fin}{Point pattern of 'case' observations classified as being inside
#' the \code{cutoff} contours of the p-value surface. An object of class
#'  \code{\link[spatstat.geom]{ppp}}.}
#' 
#' \item{fout}{Point pattern of 'case' observations falling outside the
#'  \code{cutoff} contours of the p-value surface. An object of class
#'  \code{\link[spatstat.geom]{ppp}}.}
#' 
#' \item{gin}{As \code{fin}, for the control points.}
#' 
#' \item{gout}{As \code{fout}, for the control points.}
#' 
#' \item{findex}{Numeric vector giving the raw index values of the original 
#' pattern of cases which provide \code{fin}.}
#' 
#' \item{gindex}{As \code{findex}, for the controls.}
#' 
#' \item{finsplit}{A list of the indexes in \code{findex}, with separate members
#' splitting up the indexes of case observations as falling inside each unique
#' tolerance contour.}
#' 
#' \item{ginsplit}{As \code{ginsplit}, for the controls.}
#' 
#' \item{pcmask}{The classification object of class \code{\link[spatstat.geom]{owin}}. This is
#' a pixel image mask derived from \code{pim} and \code{cutoff}.}
#' 
#' \item{pcpolys}{A list of the same length as \code{finsplit} and \code{ginsplit}, 
#' identifying each unique contour as a polygonal \code{\link[spatstat.geom]{owin}}. The order of these
#' objects in the list correspond to the membership of \code{finsplit} and \code{ginsplit}. Use
#' \code{\link{tol.classplot}} on this component to plot the classification indexes.}
#'
#'
#' @author T. M. Davies
#'
#' @references
#' 
#' Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative
#' risk, \emph{Statistics in Medicine}, \bold{29}(23) 2423-2437.
#'
#' Hazelton, M.L. and Davies, T.M. (2009), Inference based on kernel estimates
#' of the relative risk function in geographical epidemiology,
#' \emph{Biometrical Journal}, \bold{51}(1), 98-109.
#'
#' Kelsall, J.E. and Diggle, P.J. (1995), Kernel estimation of relative risk, \emph{Bernoulli},
#' \bold{1}, 3-16.
#'
#' @examples
#' 
#' data(pbc)
#' pbccas <- split(pbc)$case
#' pbccon <- split(pbc)$control
#' h0 <- OS(pbc,nstar="geometric")
#' 
#' pbcrr <- risk(pbccas,pbccon,h0=h0,tolerate=TRUE)
#' pbcclass <- tol.classify(pbcrr)
#' 
#' \dontrun{
#' plot(pbcrr)
#' points(pbcclass$fin,col="red",pch=3,cex=0.5)
#' points(pbcclass$fout,col="seagreen4",cex=0.5)
#' 
#' chrr <- risk(chorley,h0=0.7,tolerate=TRUE)
#' chclass <- tol.classify(chrr,cutoff=0.4)
#' plot(chrr,tol.args=list(levels=0.4))
#' for(i in 1:length(chclass$finsplit)){
#'    points(chrr$f$pp[chclass$finsplit[[i]]],col=i,pch=19)
#' }
#' }
#' 
#' 
#' 
#' 
#' 
#' @export
tol.classify <- function(rs,cutoff=0.05,pim=NULL,...){
  if(!inherits(rs,"rrs")) stop("'rs' argument must be of class \"rrs\"")
    
  pim <- rs$P
  fdata <- coords(rs$f$pp)
  gdata <- coords(rs$g$pp)
  
  if(is.null(pim)){
    pim <- tolerance(rs,...)
  } else{
    if(!is.im(pim)) stop("'pim' must be a pixel image of class \"im\"")
  }
  
  if(!compatible.im(rs$rr,pim)) stop("'pim' and 'rs' have incompatible pixel images -- check resolution and data")
  
  pimcut <- pim <= cutoff
  if(sum(pimcut)==0) warning("no pixel values <= 'cutoff' in p-value surface")
  
  pcmask <- owin(mask=as.matrix(pimcut),xy=list(x=pimcut$xcol,y=pimcut$yrow))
  pcpoly <- as.polygonal(pcmask)
  
  findex <- which(inside.owin(x=fdata$x,y=fdata$y,w=pcmask))
  gindex <- which(inside.owin(x=gdata$x,y=gdata$y,w=pcmask))
 
  if(length(findex)/nrow(fdata)>0.95) warning("over 95% of 'case' points inside tolerance contours")
  if(length(gindex)/nrow(gdata)>0.95) warning("over 95% of 'control' points inside tolerance contours")
  
  pareas <- summary(pcpoly)$areas
  paneg <- which(pareas<0)
  npneg <- length(paneg)
  
  if(npneg==0){
    pcpolysplit <- lapply(pcpoly$bdry,function(wn) owin(poly=list(x=wn$x,y=wn$y)))
  } else {
    papos <-  which(pareas>0)
    pholeown <- rep(NA,npneg)
    for(j in 1:npneg) pholeown[j] <- max(which(papos<paneg[j]))
    pcpolysplit <- list()
    for(i in 1:length(papos)){
      holy <- which(pholeown==i)
      temppoly <- owin(poly=pcpoly$bdry[c(papos[i],paneg[holy])])
      pcpolysplit[[i]] <- temppoly
    }
  }
  
  pclen <- length(pcpolysplit)
  fpsplit <- gpsplit <- list()
  for(i in 1:pclen){
    fpsplit[[i]] <- which(inside.owin(x=fdata$x,y=fdata$y,w=pcpolysplit[[i]]))
    gpsplit[[i]] <- which(inside.owin(x=gdata$x,y=gdata$y,w=pcpolysplit[[i]]))
  }
  
  return(list(fin=rs$f$pp[findex],
              fout=rs$f$pp[-findex],
              gin=rs$g$pp[gindex],
              gout=rs$g$pp[-gindex],
              findex=findex,
              gindex=gindex,
              finsplit=fpsplit,
              ginsplit=gpsplit,
              pcmask=pcmask,
              pcpolys=pcpolysplit))
}