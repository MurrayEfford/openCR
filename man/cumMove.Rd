\name{cumMove}
\alias{cumMove}
\alias{proportionInPolygon}

\title{

Probability Distribution After Movement

}

\description{

Compute the compounding effect of a random walk defined by a discrete kernel. The number of steps and the edge algorithm are specified by the user. The function was used to generate Fig. 3 of Efford (2022). The final distribution may be summed for points lying within an arbitrary polygon. This is a simple way to compute the expected proportion remaining within a particular region (i.e. not ``emigrating").

}

\usage{

cumMove(X, mask, kernel, edgemethod = c("truncate", "wrap", "none"), nstep = 1,
 mqarray = NULL, settlecov = NULL)

proportionInPolygon(mask, poly, cov = "pm")

}

\arguments{
  \item{X}{initial location(s) (see Details)}
  \item{mask}{habitat mask}
  \item{kernel}{kernel object}
  \item{edgemethod}{character}
  \item{nstep}{non-negative integer}
  \item{mqarray}{integer array of lookup indices}
  \item{settlecov}{character name of covariate of \code{mask}}
  \item{poly}{a polygon (see Details)}
  \item{cov}{character name of covariate of \code{mask}}
}

\details{

The input \code{X} may be -
\itemize{
\item a vector of length 2 for the coordinates of a single point 
\item a mask with covariate 'pm' representing the initial distribution
\item a SpatialPolygons object from \pkg{sp}. Animals are assumed initially to be distributed uniformly across mask points that lie within the polygon.
}

The default edgemethod truncates the kernel at the edge and re-normalizes the cell probabilities so that all destinations lie within the boundary of the mask.

\code{settlecov} may name a covariate of \code{mask} that has settlement weights in range 0--1.

For \code{proportionInPolygon}, the input mask may be the output from \code{cumMove}. The polygon \code{poly} may be specified as for \code{\link{pointsInPolygon}} (e.g., SpatialPolygons object or 2-column matrix of coordinates) or as a list with components x and y. A list of polygon specifications is also accepted.

\code{mqarray} is computed automatically if not provided. Precomputing the array can save time but is undocumented.

}

\value{

For cumMove - a mask object with initial probability distribution in covariate 'pm0' and final distribution in covariate 'pm'.

For proportionInPolygon - vector of the summed weights (probabilities) for cells centred in the polygon(s) as a proportion of all non-missing weights.

}

\seealso{

\code{\link{make.kernel}}, \code{\link{pointsInPolygon}}

}

\references{

Efford, M. G. (2022) . Efficient discretization of movement kernels for 
spatiotemporal capture--recapture. \emph{Journal of Agricultural, Biological and 
Environmental Statistics}. In press. https://doi.org/10.1007/s13253-022-00503-4

}

}
\examples{

sp <- 10
msk <- make.mask(nx = 51, ny = 51, type = 'rect', spacing = sp, 
    buffer = 0)
k <- make.kernel('BVN', 20, spacing = sp, move.a = 50, clip = TRUE, 
    sparse = TRUE)

# initial distribution a central point
X <- apply(msk, 2, mean)   
par(mfrow = c(1,4), mar = c(1,1,2,1))
for (step in 0:2) {
    X <- cumMove(X, msk, k, nstep = min(step,1))
    plot(X, cov = 'pm', dots = FALSE, legend = FALSE, breaks = 
        seq(0,0.006,0.0001))
        mtext(side = 3, line = 0, paste('Step', step), cex = 0.9)
    contour(
        x = unique(X$x), 
        y = unique(X$y), 
        z = matrix(covariates(X)$pm, nrow = length(unique(X$x))), 
        levels = c(0.0002), 
        drawlabels = FALSE,
        add = TRUE)
}

\dontrun{
# initial distribution across a polygon
X0 <- matrix(c(200,200,300,300,200,200,300,300,200,200), ncol = 2)
X <- X0
par(mfrow = c(1,4), mar = c(1,1,2,1))
for (step in 0:3) {
    X <- cumMove(X, msk, k, nstep = min(step,1))
    plot(X, cov = 'pm', dots = FALSE, legend = FALSE, breaks = 
        seq(0,0.006,0.0001))
        mtext(side = 3, line = 0, paste('Step', step), cex = 0.9)
    contour(
        x = unique(X$x), 
        y = unique(X$y), 
        z = matrix(covariates(X)$pm, nrow = length(unique(X$x))), 
        levels = c(0.0002), 
        drawlabels = FALSE,
        add = TRUE)
}
polygon(X0)
proportionInPolygon(X, X0)
}

}
