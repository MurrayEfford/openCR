\name{age.matrix}
\alias{age.matrix}

\title{Session-specific Ages}

\description{

A matrix showing the age of each animal at each secondary session (occasion).

}

\usage{

age.matrix(capthist, initialage = 0, minimumage = 0, maximumage = 1, 
    collapse = FALSE, unborn = minimumage)

}

\arguments{
  \item{capthist}{single-session capthist object}
  \item{initialage}{numeric or character name of covariate with age at first detection (optional)}
  \item{minimumage}{integer minimum age}
  \item{maximumage}{integer maximum age}
  \item{collapse}{logical; if TRUE then values for each individual are collapsed 
      as a string with no spaces}
  \item{unborn}{numeric code for age<0 }
}

\details{

\code{age.matrix} is used by \code{\link{openCR.design}} for the predictors `age' and `Age'.

Computations use the intervals attribute of \code{capthist}, which may be non-integer.

Ages are inferred for occasions before first detection, back to the minimum age.

}

\value{

Either a numeric matrix with dimensions (number of animals, number of secondary occasions)
or if \code{collapse = TRUE} a character matrix with one column.

}

\seealso{

\code{\link{openCR.design}}

}

\examples{

age.matrix(join(ovenCH), maximumage = 2, collapse = TRUE)

}

\keyword{ manip }