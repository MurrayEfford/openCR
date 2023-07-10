\name{pkernel}
\alias{gkernel}
\alias{pkernel}
\alias{dkernel}
\alias{qkernel}

\title{Kernel Distribution Functions}

\description{

Distribution of distance moved for each of the main movement kernels. Theoretical probability density, cumulative distribution function, and quantile function (inverse of the cumulative distribution function). 

}

\usage{

pkernel(q, movementmodel = c("BVN", "BVE", "BVC", "BVT", "RDE", "RDG", "RDL"), 
    move.a, move.b, truncate = Inf, lower.tail = TRUE)

dkernel(r, movementmodel = c("BVN", "BVE", "BVC", "BVT", "RDE", "RDG", "RDL"), 
    move.a, move.b, truncate = Inf)

qkernel(p, movementmodel = c("BVN", "BVE", "BVC", "BVT", "RDE", "RDG", "RDL"), 
    move.a, move.b, truncate = Inf, lower.tail = TRUE)

gkernel(r, movementmodel = c("BVN", "BVE", "BVC", "BVT", "RDE", "RDG", "RDL"),
    move.a, move.b, truncate = Inf)
    
}

\arguments{
  \item{p}{numeric vector of cumulative probabilities (0.5 for median)}
  \item{r}{numeric vector of distance moved}
  \item{q}{numeric vector of quantiles (distance moved)}
  \item{movementmodel}{character (see \link{Movement models} and \href{../doc/openCR-vignette.pdf}{openCR-vignette.pdf})}
  \item{move.a}{numeric parameter of movement kernel}
  \item{move.b}{numeric parameter of movement kernel}
  \item{truncate}{numeric q value at which distribution truncated}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].}
}

\details{

Some formulae are given in openCR-vignette.pdf. \code{gkernel} gives the 2-D probability density of the bivariate kernel \eqn{g(r) = f(r) / (2\pi r)}; the remaining functions describe the distribution of distance moved \eqn{f(r)}.

Computation of \code{qkernel} for \code{movementmodel = 'BVE'} uses 
numerical root finding (function \code{\link{uniroot}}).

Truncation (\code{truncate = limit} for finite \code{limit}) adjusts probabilities 
upwards by 1/pkernel(limit,..., truncate = Inf) so that 
pkernel(limit, ..., truncate = limit) equals 1.0.
By default the distribution is not truncated.

}

\value{

For \code{pkernel} --

Vector of cumulative probabilities corresponding to q. The cumulative probability is 1.0 for q > truncate.

For \code{dkernel} --

Vector of probability density at radial distance r (zero for r > truncate).

For \code{qkernel} --

Vector of quantiles (distances moved) corresponding to cumulative probabilities p.

For \code{gkernel} --

Vector of 2-D probability density at radial distance r (zero for r > truncate).

}

\seealso{

\link{Movement models}, \code{\link{make.kernel}}, \code{\link{matchscale}}

}

\references{

Efford, M. G. and Schofield, M. R. (2022) A review of movement models in open population capture--recapture. \emph{Methods in Ecology and Evolution} \bold{13}, 2106--2118. https://doi.org/10.1111/2041-210X.13947

}

\examples{
# plot 3 distributions chosen with matchscale to intersect at p = 0.5
q <- 0:100
plot(q, pkernel(q, 'BVN', 34), type = 'l', ylab = 'Cumulative probability')
lines(q, pkernel(q, 'BVT', move.a = 104, move.b = 5), col = 'darkgreen', lwd = 2)
lines(q, pkernel(q, 'BVT', move.a = 40, move.b = 1), col = 'orange', lwd = 2)
points(40, 0.5, pch = 16)
legend(62, 0.36, lty=1, lwd = 2, col = c('black','darkgreen','orange'), 
   legend = c('BVN sigma=34', 'BVT a=104, b=5', 'BVT a=40, b=1'))

# median
abline(v = qkernel(0.5, 'BVN', 34))
}
