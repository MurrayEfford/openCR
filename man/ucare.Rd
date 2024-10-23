\name{ucare.cjs}
\alias{ucare.cjs}

\title{Goodness-of-fit tests for the Cormack-Jolly-Seber model}

\description{

The package \pkg{R2ucare} (Gimenez et al. 2018, 2022) provides the standard tests for CJS models from Burnham et al. (1987) along with tests for multi-state models as described by Pradel et al. (2005). This function is a wrapper for the tests relevant to \pkg{openCR} (see Details). Original papers and the vignette for \pkg{R2ucare} should be consulted for interpretation.

}

\usage{

ucare.cjs(CH, tests = "all", by = NULL, verbose = TRUE, rounding = 3, ...)

}

\arguments{

  \item{CH}{capthist object suitable for openCR}
  \item{tests}{character vector with the names of specific tests (see Details) or `all'}
  \item{by}{character name of covariate in CH used to split rows of CH into separate groups }
  \item{verbose}{logical; if TRUE then additional details are tabulated}
  \item{rounding}{integer number of decimal places in output}
  \item{\dots}{other arguments passed to \code{\link[secr]{split.capthist}} if needed }
  
}

\details{

The possible tests are ``test3sr", ``test3sm", ``test2ct", ``test2cl", and ``overall_CJS".

If \code{CH} is a multi-session object then it will first be collapsed to a single-session object with \code{join} as usual in \pkg{openCR}. If \code{CH} has an intervals attribute indicating that the data are from a robust design (some intervals zero) then it will first be collapsed to one secondary session per primary session, with a warning.

If \code{by} is specified it should point to a categorical variable (factor or character) in the covariates attribute of \code{CH}. Separate tests will be conducted for each group.

\pkg{R2ucare} was removed from CRAN in May 2022, but will return at some point. In the meantime, it may be necessary to install from GitHub with

\code{
if(!require(devtools)) install.packages("devtools")
devtools::install_github("oliviergimenez/R2ucare")
}

}

\value{

A list of results, possibly nested by the grouping variable \code{by}. The verbose form includes both the overall result of each test and its breakdown into components (`details').

}

\references{

Burnham, K. P., Anderson, D. R., White, G. C., Brownie, C. and Pollock, K. H. (1987) \emph{Design and Analysis Methods for Fish Survival Experiments Based on Release-Recapture}. American Fisheries Society Monograph 5. Bethesda, Maryland, USA.

Choquet, R., Lebreton, J.-D., Gimenez, O., Reboulet, A.-M. and Pradel, R. (2009) U-CARE: Utilities for performing goodness of fit tests and manipulating CApture-REcapture data. \emph{Ecography} \bold{32}, 1071--1074.

Gimenez, O., Lebreton, J.-D., Choquet, R. and Pradel, R. (2018) R2ucare: An R package to perform 
goodness-of-fit tests for capture--recapture models. \emph{Methods in Ecology and Evolution} \bold{9}, 1749--1754.

Gimenez, O., Lebreton, J.-D., Choquet, R. and Pradel, R. (2022) R2ucare: Goodness-of-Fit Tests
  for Capture-Recapture Models. R package version 1.0.2. \url{https://github.com/oliviergimenez/R2ucare/}

Lebreton, J.-D., Burnham, K. P., Clobert, J., and Anderson, D. R. (1992) Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. \emph{Ecological Monographs} \bold{62}, 67--118.

Pradel, R., Gimenez O. and Lebreton, J.-D. (2005) Principles and interest of GOF tests for multistate capture--recapture models. \emph{Animal Biodiversity and Conservation} \bold{28}, 189--204.

}

\seealso{

\code{\link{m.array}}

}

\examples{

if (requireNamespace("R2ucare"))
    ucare.cjs(dipperCH, verbose = FALSE, by = 'sex')

}

\keyword{ htest }
