\name{Movement models}
\alias{Movement models}

\title{ List of Movement Models }

\description{

Movement of activity centres between primary sessions is modelled in \pkg{openCR} as a random walk with step length governed by a circular probability kernel. The argument `movementmodel' defines the kernel in several functions. More detail is provided in the vignettes \href{../doc/openCR-vignette.pdf}{openCR-vignette.pdf}.
}

\section{Movement models in \pkg{openCR} 2.2}{

Kernel models:

 \tabular{lllll}{
 Kernel \tab  \tab Description \tab  \tab Parameters \cr
 BVN    \tab  \tab bivariate normal \tab  \tab move.a\cr
 BVE    \tab  \tab bivariate Laplace \tab  \tab move.a\cr
 BVC    \tab  \tab bivariate Cauchy distribution \tab  \tab move.a \cr
 BVT    \tab  \tab bivariate t-distribution (2Dt of Clark et al. 1999) \tab  \tab move.a, move.b \cr
 RDE    \tab  \tab exponential distribution of distance moved cf Ergon and Gardner (2014) \tab  \tab move.a\cr
 RDG    \tab  \tab gamma distribution of distance moved cf Ergon and Gardner (2014) \tab  \tab move.a, move,b\cr
 RDL    \tab  \tab log-normal distribution of distance moved cf Ergon and Gardner (2014) \tab  \tab  move.a, move.b\cr
 RDLS*   \tab\tab log-sech distribution of distance moved (Van Houtan et al. 2007) \tab\tab move.a, move.b \cr
 UNI    \tab  \tab uniform within kernel radius, zero outside \tab  \tab (none) \cr
 BVNzi  \tab  \tab zero-inflated BVN \tab  \tab move.a, move.b\cr
 BVEzi  \tab  \tab zero-inflated BVE \tab  \tab move.a, move.b\cr
 RDEzi  \tab  \tab zero-inflated RDE \tab  \tab move.a, move.b\cr
 UNIzi  \tab  \tab zero-inflated UNI \tab  \tab move.a \cr
 }
 * incomplete implementation

Kernel-free models (buffer dependent):
 \tabular{lllll}{
 Model \tab  \tab Description \tab  \tab Parameters \cr
 IND \tab  \tab independent relocation within habitat mask (Gardner et al. 2018) \tab  \tab (none) \cr
 INDzi \tab  \tab zero-inflated IND \tab  \tab move.a \cr
 }

}

\section{Relationships among models}{

Some models may be derived as special cases of others, for example

\tabular{lllll}{
 General           \tab\tab Condition \tab\tab Equivalent to \cr
 BVT \tab\tab large move.b (df \eqn{\infty}) \tab\tab BVN \cr
 BVT \tab\tab move.b = 0.5 (df 1)  \tab\tab BVC \cr 
 RDG \tab\tab move.b = 1    \tab\tab RDE \cr 
 RDG \tab\tab move.b = 2    \tab\tab BVE \cr
 BVNzi \tab\tab large move.a \tab\tab UNIzi \cr
}

RDL and RDG are almost indistinguishable when move.b > 2.

}

\section{Deprecated names of movement models}{

These old names appeared in earlier releases. They still work, but may be removed in future.
 \tabular{lll}{
 Old \tab \tab New\cr
 normal \tab  \tab BVN\cr
 exponential \tab  \tab BVE\cr
 t2D \tab  \tab BVT \cr
 frE \tab  \tab RDE \cr
 frG \tab  \tab RDG \cr
 frL \tab  \tab RDL \cr
 uniform \tab  \tab UNI \cr
 frEzi \tab  \tab RDEzi\cr
 uniformzi \tab  \tab UNIzi \cr
 }

}

\section{Additional movement models that may be removed without notice}{

\tabular{lllll}{
Kernel \tab \tab Description \tab \tab Parameters \cr
annular \tab \tab non-zero only at centre and edge cells (after clipping at kernelradius) \tab \tab move.a \cr
annularR \tab \tab non-zero only at centre and a ring of cells at radius R \tab \tab move.a, move.b \cr
}

``annularR'' uses a variable radius (R = move.b x kernelradius x spacing) and weights each cell according to the length of arc it intersects; ``annularR'' is not currently allowed in \code{openCR.fit}. For the `annular' models 'move.a' is the proportion at the centre (probability of not moving). 
}

\seealso{

\code{\link{make.kernel}}, 
\code{\link{gkernel}}, 
\code{\link{dkernel}}, 
\code{\link{pkernel}}, 
\code{\link{qkernel}},
\code{\link{openCR.fit}}

}

\references{

Clark, J. S, Silman, M., Kern, R., Macklin, E. and HilleRisLambers, J. (1999) Seed dispersal near and far: patterns across temperate and tropical forests. \emph{Ecology} \bold{80}, 1475--1494.

Efford, M. G. and Schofield, M. R. (2022) A review of movement models in open population capture--recapture. \emph{Methods in Ecology and Evolution} \bold{13}, 2106--2118. https://doi.org/10.1111/2041-210X.13947

Ergon, T. and Gardner, B. (2014) Separating mortality and emigration: modelling space use, dispersal and survival with robust-design spatial capture--recapture data. \emph{Methods in Ecology and Evolution} \bold{5}, 1327--1336.

Gardner, B., Sollmann, R., Kumar, N. S., Jathanna, D. and Karanth, K. U. (2018) State space and movement specification in open population spatial capture--recapture models. \emph{Ecology and Evolution} \bold{8}, 10336--10344 \doi{10.1002/ece3.4509}.

Nathan, R., Klein, E., Robledo-Arnuncio, J. J. and Revilla, E. (2012) Dispersal kernels: review. In: J. Clobert et al. (eds) \emph{Dispersal Ecology and Evolution}. Oxford University Press. Pp. 187--210.

Van Houtan, K. S., Pimm, S. L., Halley, J. M., Bierregaard, R. O. Jr and Lovejoy, T. E. (2007) Dispersal of Amazonian birds in continuous and fragmented forest. \emph{Ecology Letters} \bold{10,} 219--229.

}