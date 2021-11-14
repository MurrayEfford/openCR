\name{classMembership}
\alias{classMembership}
\alias{classMembership.openCR}

\title{Class Membership Probability for Mixture Models}

\description{

Finite mixture models treat class membership as a latent random variable. The probability of an individual's membership in each class may be inferred retrospectively from the relative likelihoods.

}

\usage{

\method{classMembership}{openCR}(object, fullCH = NULL, ...)

}

\arguments{

  \item{object}{fitted model of class openCR}
  \item{fullCH}{capthist object (optional)}
  \item{\dots}{other arguments (not used)}
  
}

\details{

It is assumed that the input model includes finite mixture terms (h2 or h3).

As the detection histories are saved in compressed (``squeezed'') form in openCR objects the original animal identifiers are lost and the order of animals may change. These may be restored by providing \code{fullCH}.

No class can be assigned from a CJS model for animals detected only in the final session.

}

\value{

Matrix with one row per individual and columns for each class and the class number of the most likely class.

}

\seealso{

\code{\link{openCR.fit}}, 
\code{\link{squeeze}}

}

\note{

In earlier versions \code{\link{openCR.fit}} always computed class membership and saved it in component `posterior' of the fitted model. \code{classMembership} replaces that functionality.

}

\examples{

\dontrun{
jch <- join(ovenCH)   
fit <- openCR.fit(ovenCH, model=p~h2)
classMembership(fit, jch)
}

}

\keyword{models}
