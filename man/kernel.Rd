\name{make.kernel}
\alias{make.kernel}
\alias{plot.kernel}
\alias{summary.kernel}

\title{Discrete Movement Kernel}

\description{

Functions to create, plot and summarise a discrete representation of a movement kernel.

}

\usage{

make.kernel(movementmodel = c("BVN", "BVE", "BVC", "BVT","RDE", "RDG", "RDL", "UNI"), 
    kernelradius = 10, spacing, move.a, move.b, 
    sparsekernel = FALSE, clip = FALSE, normalize = TRUE,     
    stat = c('estimate','lcl', 'ucl'), session = 1, r0 = 1/sqrt(pi), ...)

\method{plot}{kernel}(x, type = "kernel", contour = FALSE, levels = NULL, text = FALSE,
    title = NULL, add = FALSE, xscale = 1, ...)

\method{summary}{kernel}(object, ...)

}
  
\arguments{

  \item{movementmodel}{character or function or openCR object}
  \item{kernelradius}{integer radius of kernel in grid cells}
  \item{spacing}{numeric spacing between cell centres}
  \item{move.a}{numeric parameter of kernel}
  \item{move.b}{numeric parameter of kernel}
  \item{sparsekernel}{logical; if TRUE then only cardinal and intercardinal axes are included}
  \item{clip}{logical; if TRUE then corner cells are removed}
  \item{normalize}{logical; if TRUE then cell values are divided by their sum}
  \item{stat}{character; predicted statistic to use for move.a (openCR object only) }
  \item{session}{integer; session for move.a, move.b if input is fitted model}
  \item{r0}{numeric; effective radius of zero cell for movement models}
  
  \item{x}{kernel object from \code{make.kernel}}
  \item{type}{character; plot style (see Details)}
  \item{contour}{logical; if TRUE then contour lines are overlaid on any plot}
  \item{levels}{numeric vector of contour levels}
  \item{text}{logical; if TRUE then cell probabilities are overprinted, rounded to 3 d.p.}
  \item{title}{character; if NULL a title is constructed automatically}
  \item{add}{logical; if TRUE a line is added to an existing plot (types "gr", "fr", "Fr")}
  \item{xscale}{numeric multiplier for distance axis (0.001 for distances in km)}
  \item{\dots}{other arguments passed to \code{\link{predict.openCR}} (\code{make.kernel}) or
      \code{\link{plot.mask}} (plot type "kernel") or \code{\link{lines}} ( plot types "gr", "fr", "Fr") (not used by summary method)}

  \item{object}{kernel object from \code{make.kernel}}

}

\details{

A kernel object is a type of mask with cell probabilities stored in the covariate `kernelp'. All kernels are truncated at kernelradius x spacing.

The \code{movementmodel} may also be a function or a previously fitted openCR model that includes movement. If a fitted openCR object, parameter values and kernel attributes are derived from that object and other arguments are ignored.

The parameter `move.a' is a scale parameter in metres, except for the UNIzi and INDzi models for which it is the zero-inflation parameter (`move.b' is the zero-inflation parameter for BVNzi, BVEzi and RDEzi). 

'Sparse' kernels include only those grid cells that lie on 4 axes (N-S, E-W, NW-SE, NE-SW); cell probabilities are adjusted to maintain nearly the same distance distribution as the non-sparse equivalents.

Movement models are listed in \link{Movement models} and further described in the vignettes \href{../doc/openCR-vignette.pdf}{openCR-vignette.pdf}.

Plot type may be one or more of --
\tabular{lll}{
`kernel' \tab \tab coloured 2-D depiction \cr
`gr' \tab \tab cross-section through the origin of \eqn{g(r)} (the 2-D kernel) \cr
`fr' \tab \tab continuous probability density \eqn{f(r)} \cr
`Fr' \tab \tab cumulative probability distribution \eqn{F(r)} \cr
}

Type ``kernel" by default includes an informative title with font size from the graphical parameter `cex.main'. Set \code{title = ""} to suppress the title.

Useful properties of theoretical (not discretized) kernels may be recovered with \code{\link{matchscale}}, \code{\link{pkernel}}, \code{\link{dkernel}} and \code{\link{qkernel}}.

The obscure argument \code{r0} controls the value assigned to the central cell of a discretized kernel. For positive \code{r0} the value is F(r0*cellsize), where F is the cumulative probability distribution of distance moved. Otherwise the cell is assigned the value g(0)*cellarea, where g() is the 2-D kernel probability density (this fails where g(0) is undefined or infinite).

}

\value{

\code{make.kernel} returns an object of class c('kernel','mask','data.frame').

The kernel object has attributes: 

\tabular{ll}{
Attribute \tab Description \cr
movementmodel \tab saved input \cr
K2 \tab saved kernelradius \cr
move.a \tab saved input \cr
move.b \tab saved input \cr
distribution \tab empirical cumulative distribution function \cr
}

The empirical cumulative distribution is a dataframe with columns for the sorted cell radii `r' and the 
associated cumulative probability `cumprob' (one row per cell).

\code{summary.kernel} returns an object with these components, displayed with the corresponding print method.

\tabular{ll}{
Component \tab Description \cr
k2        \tab kernel radius in mask cells \cr
spacing   \tab cell width \cr
ncells    \tab number of cells in kernel \cr
movementmodel \tab movement model code \cr
move.a    \tab first (scale) parameter \cr
move.b    \tab second (shape) parameter \cr
mu        \tab mean of logs (RDL only; from move.a) \cr
s         \tab SD of logs (RDL only; from move.b) \cr
expectedmove \tab mean movement (untruncated) \cr
expectedmovetr \tab mean movement (trucated at kernel radius) \cr
expectedmoveemp \tab mean computed directly from kernel cell values as sum(r.p) \cr
ptruncated    \tab proportion of theoretical distribution truncated at radius \cr
expectedq50   \tab theoretical (untruncated) median \cr
expectedq90   \tab theoretical (untruncated) 90th percentile \cr
expectedq50tr \tab theoretical truncated median \cr
expectedq90tr \tab theoretical truncated 90th percentile\cr
}

The empirical mean in \code{expectedmoveemp} is usually the most pertinent property of a fitted kernel.

}


\note{

The plot method for kernels supercedes the function \code{plotKernel} that has been removed.

}

\references{

Clark, J. S, Silman, M., Kern, R., Macklin, E. and HilleRisLambers, J. (1999) Seed dispersal near and far: patterns across temperate and tropical forests. \emph{Ecology} \bold{80}, 1475--1494.

Efford, M. G. and Schofield, M. R. (2022) A review of movement models in open population 
capture--recapture. \emph{Methods in Ecology and Evolution}. In press. 
https://doi.org/10.1111/2041-210X.13947

Ergon, T. and Gardner, B. (2014) Separating mortality and emigration: modelling space use, dispersal and survival with robust-design spatial capture--recapture data. \emph{Methods in Ecology and Evolution} \bold{5}, 1327--1336.

Nathan, R., Klein, E., Robledo-Arnuncio, J. J. and Revilla, E. (2012) Dispersal kernels: review. In: J. Clobert et al. (eds) \emph{Dispersal Ecology and Evolution}. Oxford University Press. Pp. 187--210.
 
}

\seealso{

\link{Movement models}, \code{\link{mask}}, \code{\link{matchscale}}, \code{\link{dkernel}}, \code{\link{pkernel}}, \code{\link{qkernel}}

}

\examples{

k <- make.kernel(movementmodel = 'BVT', spacing = 10, move.a = 20, move.b = 1)
summary(k)

# read a previously fitted movement model packaged with 'openCR'
fit <- readRDS(system.file("exampledata", "spmOV.RDS", package = "openCR"))
k <- make.kernel(fit)
plot(k)
if (interactive()) {
   spotHeight(k, dec = 3)  # click on points; Esc to exit
}
}

\keyword{ hplot }
