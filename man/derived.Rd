\name{derived}
\alias{derived.openCR}
\alias{derived.openCRlist}
\alias{openCR.esa}
\alias{openCR.pdot}

\title{

Derived Parameters From openCR Models

}

\description{

For ..CL openCR models, compute the superpopulation size or density. 
For all openCR models, compute the
time-specific population size or density from the estimated
superpopulation size and the turnover parameters.

}

\usage{

\method{derived}{openCR}(object, newdata = NULL, all.levels = FALSE, Dscale = 1, 
    HTbysession = FALSE, ...)
\method{derived}{openCRlist}(object, newdata = NULL, all.levels = FALSE, Dscale = 1, 
    HTbysession = FALSE, ...)
openCR.esa(object, bysession = FALSE, stratum = 1)
openCR.pdot(object, bysession = FALSE, stratum = 1)

}

\arguments{
  \item{object}{fitted openCR model}
  \item{newdata}{optional dataframe of values at which to evaluate model}
  \item{all.levels}{logical; passed to \code{\link{makeNewData}} if newdata not specified}
 % \item{all.levels}{logical; passed to \code{\link{openCR.make.newdata}} if newdata not specified}
  \item{Dscale}{numeric to scale density}
  \item{HTbysession}{logical; Horvitz-Thompson estimates by session (see Details)}
  \item{...}{other arguments (not used)}
  \item{bysession}{logical; if TRUE then esa or pdot is computed separately for each session}
  \item{stratum}{integer}
}

\details{

Derived estimates of density and superD are multiplied by \code{Dscale}. Use \code{Dscale = 1e4} for animals per 100 sq. km.
\code{openCR.esa} and \code{openCR.pdot} are used internally by \code{derived.openCR}.

If \code{HTbysession} then a separate H-T estimate is derived for each primary session; otherwise a H-T estimate of the superpopulation is used in combination with turnover parameters (phi, beta) to obtain session-specific estimates. Results are often identical.

The output is an object with its own print method (see \code{\link{print.derivedopenCR}}).

The code does not yet allow user-specified newdata.

}

\value{

  \code{derived} returns an object of class c(``derivedopenCR",``list"), list with these components:
  \item{totalobserved}{number of different individuals detected}
  \item{parameters}{character vector; names of parameters in model (excludes derived parameters)}
  \item{superN }{superpopulation size (non-spatial models only)}
  \item{superD}{superpopulation density (spatial models only)}
  \item{estimates}{data frame of counts and estimates}
  \item{Dscale}{numeric multiplier for printing densities}

If \code{newdata} has multiple levels then the value is a list of such objects, one for each level.

  \code{openCR.pdot} returns a vector of experiment-wide detection
  probabilities under the fitted model (one for each detected animal).

  \code{openCR.esa} returns a vector of effective sampling areas under
  the fitted model (one for each detected animal). If 'bysession = TRUE' the 
  result is a list with one component per session.
  
}

\seealso{

  \code{\link{openCR.fit}}, \code{\link{print.derivedopenCR}}
  
}

\note{

Prior to 1.4.5, openCR.esa did not expand the result for squeezed capture 
histories (freq>1) and did not return a list when bysession = TRUE.

}

\examples{

\dontrun{

# override default method to get true ML for L1
L1CL <- openCR.fit(ovenCH, type = 'JSSAlCL', method = 'Nelder-Mead')
predict(L1CL)
derived(L1CL)

## compare to above
L1 <- openCR.fit(ovenCH, type = 'JSSAl', method = 'Nelder-Mead')
predict(L1)
derived(L1)

}

}

\keyword{ models }
