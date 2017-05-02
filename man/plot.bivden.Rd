\name{plot.bivden}
\Rdversion{1.1}
\alias{plot.bivden}
\alias{plot.rrs}
\alias{plot.msden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting a bivariate kernel density/intensity, relative risk surface, or multi-scale density/intensity
}
\description{
\code{plot} methods for classes \code{"\link{bivden}"}, \code{"\link{rrs}"}, and \code{"\link{msden}"}
}
\usage{
\method{plot}{bivden}(x, what = c("z", "edge", "bw"), add.pts = FALSE, auto.axes = TRUE, \dots)

\method{plot}{rrs}(x, auto.axes = TRUE, \dots)

\method{plot}{msden}(x, what = c("z", "edge", "bw"), sleep = 0.2, override.par = TRUE, \dots)
}
\arguments{
  \item{x}{
An object of class \code{\link{bivden}}, \code{\link{rrs}}, or \code{\link{msden}}.
}
\item{what}{
A character string to select plotting of result (\code{"z"}; default); edge-correction surface (\code{"edge"}); or variable bandwidth surface (\code{"bw"}).
}
\item{add.pts}{
Logical value indicating whether to add the observations to the image plot using default \code{\link{points}}.
}
\item{auto.axes}{
Logical value indicating whether to display the plot with automatically added x-y axes and an `L' box in default styles.
}
\item{sleep}{
Single positive numeric value giving the amount of time (in seconds) to \code{\link[base]{Sys.sleep}} before drawing the next image in the animation.
}
\item{override.par}{
Logical value indicating whether to override the existing graphics device parameters prior to plotting, resetting \code{mfrow} and \code{mar}. See `Details' for when you might want to disable this.
}
  \item{...}{
Additional graphical parameters to be passed to \code{\link[spatstat]{plot.im}}, or in one instance, to \code{\link[spatstat]{plot.ppp}} (see `Details').
}
}

  
\details{
In all instances, visualisation is deferred to \code{\link[spatstat]{plot.im}}, for which there are a variety of customisations available the user can access via \code{...}. The one exception is when plotting observation-specific \code{"diggle"} edge correction factors---in this instance, a plot of the spatial observations is returned with size proportional to the influence of each correction weight.

When plotting a \code{\link{rrs}} object, a pre-computed \emph{p}-value surface (see argument \code{tolerate} in \code{\link{risk}}) will automatically be superimposed at a significance level of 0.05. Greater flexibility in visualisation is gained by using \code{\link{tolerance}} in conjunction with \code{\link{contour}}.

The \code{\link{msden}} object is plotted as an animation, one pixel image after another, separated by \code{sleep} seconds. If instead you intend the individual images to be plotted in an array of images, you should first set up your plot device layout, and ensure \code{override.par = FALSE} so that the function does not reset these device parameters itself. In such an instance, one might also want to set \code{sleep = 0}.
}

\value{
Plots to the relevant graphics device.
}

%\references{
%}

\author{
T.M. Davies
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# to be filled 
}

