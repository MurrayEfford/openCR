\name{expected.d}
\alias{expected.d}

\title{Expected Distance Moved}

\description{

Movement models in \pkg{openCR} differ in their parameterisation so direct 
comparison can be difficult. The expected distance moved is a convenient 
statistic common to all models. This function computes the expected distance 
from various inputs, including fitted models.

}

\usage{

expected.d(movementmodel, move.a, move.b, truncate = Inf, mask = NULL, 
    min.d = 1e-4, ...)

}

\arguments{

  \item{movementmodel}{character or function or kernel or openCR object}
  \item{move.a}{numeric parameter of kernel}
  \item{move.b}{numeric parameter of kernel}
  \item{truncate}{radius of truncation}
  \item{mask}{habitat mask object}
  \item{min.d}{numeric lower bound of integration (see Details)}
  \item{\dots}{other arguments passed to \code{\link{make.kernel}} if input is a fitted model}

}

\details{

The input \code{movementmodel} may be
\itemize{
\item fitted openCR model
\item user kernel function g(r)
\item kernel object
\item character name of kernel model see \link{Movement models} 
}

If \code{truncate} (R) is finite or \code{movementmodel} is a function then the 
expected value is computed by numerical integration \eqn{E(d) = \int_0^R r.f(r) dr}. 
In the event that f(0) is not finite, \code{min.d} is used as the lower bound. 

\code{mask} is used only for `uncorrelated' and `uncorrelatedzi' movement. 
For these models the expected movement is merely the average distance 
between points on the mask, weighted by (1-zi) if zero-inflated (uncorrelatedzi).
 
The \dots argument is useful for (i) selecting a session from a fitted model, or
(ii) specifying the upper or lower confidence limits from a single-parameter 
fitted model via the `stat' argument of \code{\link{make.kernel}}.

}

\value{

A numeric value (zero for 'static' model, NA if model unrecognised).

}

\seealso{

\link{Movement models}, \code{\link{make.kernel}}, \code{\link{pkernel}}, \code{\link{qkernel}}

}

\examples{

expected.d('BVT', move.a = 20, move.b = 1)
expected.d('BVT', move.a = 20, move.b = 1, truncate = 300)

k <- make.kernel(movementmodel = 'BVT', spacing = 10, move.a = 20, move.b = 1)
expected.d(k)

}
