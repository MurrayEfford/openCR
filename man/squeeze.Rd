\name{squeeze}
\alias{squeeze}
\alias{unsqueeze}

\title{Unique Capture Histories}

\description{

Compresses or expands capthist objects.

}

\usage{

squeeze(x)
unsqueeze(x)

}

\arguments{
  \item{x}{secr capthist object}
  }
  
\details{

Although \code{squeeze} may be applied to spatial capthist objects, the
effect is often minimal as most spatial histories are unique.

The `freq' covariate is used by \code{openCR.fit} to weight summaries and likelihoods. It is currently ignored by \code{secr.fit}.

}

\value{

Both functions return a capthist object with one row for each unique capture history (including covariates). The individual covariate `freq' records the number of instances of each unique history in the input.

}

\seealso{

\code{\link{openCR.fit}}

}

\examples{

squeeze(captdata)

}

\keyword{ manip }