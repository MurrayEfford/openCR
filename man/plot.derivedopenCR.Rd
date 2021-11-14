\name{plot.derivedopenCR}
\alias{plot.derivedopenCR}

\title{Plot Derived Estimates}

\description{

Session-specific estimates of the chosen parameter are plotted.

}

\usage{

\method{plot}{derivedopenCR} (x, par = "phi", add = FALSE, xoffset = 0, ylim = NULL, 
    useintervals = TRUE, intermediate.x = TRUE, ...)
    
}

\arguments{
  \item{x}{openCR object from openCR.fit}
  \item{par}{character names of parameter to plot}
  \item{add}{logical; if TRUE then points are added to an existing plot}
  \item{xoffset}{numeric offset to be added to all x values}
  \item{ylim}{numeric vector of limits on y-axis}
  \item{useintervals}{logical; if TRUE then x values are spaced according to 
        the intervals attribute}
  \item{intermediate.x}{logical; if TRUE then turnover parameters are plotted at 
        the mid point on the x axis of the interval to which they relate}
  \item{\dots}{other arguments passed to \code{\link{plot}}, \code{\link{points}} and \code{\link{segments}}}
}

\details{

If ylim is not provided it is set automatically.

Confidence intervals are not available in this version.

}

\value{

The x coordinates (including xoffset) are returned invisibly.

}

\seealso{

\code{\link{plot.openCR}}

}

\examples{

\dontrun{

fit <- openCR.fit(dipperCH, type='JSSAfCL', model = phi~session)
der <- derived(fit)
plot(der,'N', pch = 16, cex = 1.3)

}

}

\keyword{ hplot }