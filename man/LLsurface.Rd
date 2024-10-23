\name{LLsurface}
\alias{LLsurface}
\alias{LLsurface.openCR}

\title{Plot Likelihood Surface}

\description{

Calculate log likelihood over a grid of values of two beta parameters
from a fitted openCR model and optionally make an approximate contour plot of the log
likelihood surface.

This is a method for the generic function \code{LLsurface} defined in \pkg{secr}.

}

\usage{

\method{LLsurface}{openCR}(object, betapar = c("phi", "sigma"), xval = NULL, yval = NULL, 
   centre = NULL, realscale = TRUE, plot = TRUE, plotfitted = TRUE, ncores = NULL, ...)

}

\arguments{
  \item{object}{\code{openCR} object output from \code{openCR.fit}}
  \item{betapar}{character vector giving the names of two beta parameters}
  \item{xval}{vector of numeric values for x-dimension of grid}
  \item{yval}{vector of numeric values for y-dimension of grid}
  \item{centre}{  vector of central values for all beta parameters}
  \item{realscale}{logical. If TRUE input and output of x and y is on
        the untransformed (inverse-link) scale.}
  \item{plot}{logical. If TRUE a contour plot is produced}
  \item{plotfitted}{logical. If TRUE the MLE from \code{object} is shown on the plot (+)}
  \item{ncores}{integer number of cores available for parallel processing}
  \item{\dots}{other arguments passed to \code{\link{contour}}}
}

\details{

\code{centre} is set by default to the fitted values of the beta
parameters in \code{object}. This has the effect of holding parameters
other than those in \code{betapar} at their fitted values.

If \code{xval} or \code{yval} is not provided then 11 values are set at
equal spacing between 0.8 and 1.2 times the values in \code{centre} (on
the `real' scale if \code{realscale} = TRUE and on the `beta' scale
otherwise).

Contour plots may be customized by passing graphical parameters through
the \dots argument.

The value of \code{ncores} is passed to \code{openCR.fit}.

}

\value{

Invisibly returns a matrix of the log likelihood evaluated at each
grid point

}

\seealso{

\code{\link[secr]{LLsurface.secr}}

}

\note{

  \code{LLsurface.openCR} works for named `beta' parameters rather than
  `real' parameters. The default \code{realscale = TRUE} only works for
  beta parameters that share the name of the real parameter to which
  they relate i.e. the beta parameter for the base level of the real
  parameter. This is because link functions are defined for real
  parameters not beta parameters.
  
  Handling of multiple threads was changed in version 1.5.0 to align 
  with \code{\link[secr]{LLsurface.secr}}.

  The contours are approximate because they rely on interpolation. 
  
}

\examples{

# not yet

}

\keyword{ hplot }