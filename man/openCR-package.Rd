\name{openCR-package}
\alias{openCR-package}
\alias{openCR}

\docType{package}

\title{Open Population Capture--Recapture Models }

\description{
  
  Functions for non-spatial open population analysis by
  Cormack-Jolly-Seber (CJS) and Jolly-Seber-Schwarz-Arnason (JSSA)
  methods, and by spatially explicit extensions of these
  methods. The methods build on Schwarz and Arnason (1996), Borchers and 
  Efford (2008) and Pledger et al. (2010) (see \href{../doc/openCR-vignette.pdf}{vignette} 
  for more comprehensive references and likelihood). The parameterisation of JSSA 
  recruitment is flexible (options include population growth rate \eqn{\lambda}, 
  per capita recruitment \eqn{f} and seniority \eqn{\gamma}). Spatially explicit 
  analyses may assume home-range centres are fixed or allow dispersal between
  primary sessions according to various probability kernels, including bivariate 
  normal (BVN) and bivariate \emph{t} (BVT) (Efford and Schofield 2022).

}

\details{

\tabular{ll}{
Package: \tab openCR\cr
Type: \tab Package\cr
Version: \tab 2.2.6\cr
Date: \tab 2023-07-08\cr
License: \tab GNU General Public License Version 2 or later\cr
}

Data are observations of marked individuals from a `robust' sampling
design (Pollock 1982). Primary sessions may include one or more
secondary sessions. Detection histories are assumed to be stored in an object of class
`capthist' from the package \pkg{secr}. Grouping of occasions into
primary and secondary sessions is coded by the `intervals' attribute
(zero for successive secondary sessions).

A few test datasets are provided (\code{microtusCH}, \code{FebpossumCH}, \code{dipperCH}, 
\code{gonodontisCH}, \code{fieldvoleCH}) and some from \pkg{secr} are also suitable 
e.g. \code{ovenCH} and \code{OVpossumCH}.

Models are defined using symbolic formula notation. Possible predictors
include both pre-defined variables (b, session etc.), corresponding to
`behaviour' and other effects), and user-provided covariates.

Models are fitted by numerically maximizing the likelihood. The function
\code{\link{openCR.fit}} creates an object of class
\code{openCR}. Generic methods (print, AIC, etc.) are provided
for each object class.

A link at the bottom of each help page takes you to the help index. 

See \href{../doc/openCR-vignette.pdf}{openCR-vignette.pdf} for more.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\references{

  Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
  likelihood methods for capture--recapture studies. \emph{Biometrics}
  \bold{64}, 377--385.

  Efford, M. G. and Schofield, M. R. (2020) A spatial open-population capture--recapture model.
  \emph{Biometrics} \bold{76}, 392--402.
  
  Efford, M. G. and Schofield, M. R. (2022) A review of movement models in open population 
capture--recapture. \emph{Methods in Ecology and Evolution}. In press. 
https://doi.org/10.1111/2041-210X.13947

  Glennie, R., Borchers, D. L., Murchie, M. Harmsen, B. J., and Foster, R. J. (2019) Open population maximum likelihood spatial capture--recapture. \emph{Biometrics} \bold{75}, 1345--1355

  Pledger, S., Pollock, K. H. and Norris, J. L. (2010) Open
  capture--recapture models with heterogeneity: II. Jolly-Seber
  model. \emph{Biometrics} \bold{66}, 883--890.

  Pollock, K. H. (1982) A capture--recapture design robust to unequal 
  probability of capture. \emph{Journal of Wildlife Management} \bold{46}, 752--757.
  
  Schwarz, C. J. and Arnason, A. N. (1996) A general methodology for the
  analysis of capture-recapture experiments in open
  populations. \emph{Biometrics} \bold{52}, 860--873.

}

\keyword{ package }

\seealso{
  
  \code{\link{openCR.fit}}, \code{\link[secr]{capthist}}, \code{\link{ovenCH}}
  
}
\examples{

\dontrun{

## a CJS model is fitted by default
openCR.fit(ovenCH)

}

}
