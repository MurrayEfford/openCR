\name{matchscale}
\alias{matchscale}

\title{Match Kernel}

\description{

Finds scale parameter (move.a) of a movement model that corresponds to desired quantile, or expected distance moved.

}

\usage{

matchscale(movementmodel, q = 40, expected = NULL, p = 0.5, lower = 1e-05, upper = 1e+05, 
   move.b = 1, truncate = Inf)

}

\arguments{
  \item{movementmodel}{character (see \link{Movement models} and \href{../doc/openCR-vignette.pdf}{openCR-vignettes.pdf})}
  \item{q}{desired quantile (distance moved)}
  \item{expected}{numeric expected distance moved}
  \item{p}{cumulative probability}
  \item{move.b}{shape parameter of movement kernel}
  \item{lower}{lower bound interval to search }
  \item{upper}{upper bound interval to search }
  \item{truncate}{numeric q value at which distribution truncated}
}

\details{

The default behaviour is to find the movement parameter for the given combination of q and p.

The alternative, when a value is provided for `expected', is to find the movement parameter 
corresponding to the given expected distance.

The \code{truncate} argument must be specified for movementmodel `UNIzi`. For movementmodel `UNI' there is no parameter and the radius of truncation is varied to achieve the requested quantile q corresponding to cumulative probability p, or the desired expected distance.

}

\value{

Numeric value for move.a (scale parameter or zero-inflation in the case of `UNIzi') or truncation radius (`UNI').

}

\seealso{

\link{Movement models}, \code{\link{pkernel}}, \code{\link{make.kernel}}, \code{\link{expected.d}}

}

\examples{

matchscale('BVN', 40, 0.5)
matchscale('BVT', 40, 0.5, move.b = 1)
matchscale('BVT', 40, 0.5, move.b = 5)
matchscale('BVT', move.b = 5, expected = 10)

}
