\name{strata}
\alias{strata}
\alias{strata<-}

\title{ Stratum names }

\description{
  Extract or replace the stratum names of a \code{capthist} object.
}

\usage{
strata(object, \dots)
strata(object) <- value
}

\arguments{
  \item{object}{ object with `stratum' attribute e.g. \code{capthist} }
  \item{value}{ character vector or vector that may be coerced to character, one value per stratum }
  \item{\dots}{ other arguments (not used) }
}

\details{
  Replacement values will be coerced to character.
}

\value{

  a character vector with one value for each session in \code{capthist}.

}

\note{

\pkg{openCR} uses the term `stratum' for an independent set of samples, 
rather like a `session' in \pkg{secr}. Strata offer flexibility in defining 
and evaluating between-stratum models. The log likelihood for a stratified
model is the sum of the separate stratum log likelihoods. Although this 
assumes independence of sampling, parameters may be shared across strata, 
or stratum-specific parameter values may be functions of stratum-level covariates. 
The detector array and mask can be specified separately for each stratum. 

For open population analyses, each stratum comprises both primary and secondary 
sessions of Pollock's robust design `joined' in a single-session capthist object. 

The function \code{\link{stratify}} can be useful for manipulating data into 
multi-stratum form.

Models are stratified only if the argument \code{stratified} of 
\code{openCR.fit()} is set to TRUE. Strata will otherwise be treated as 
primary sessions and concatenated as usual with \code{join()}.

}

\seealso{ 

  \code{\link{openCR.fit}}, \code{\link{session}}, \code{\link{stratify}}

}

\examples{
  # artificial example, treating years as strata
  strata(ovenCH)
}

\keyword{ models }

