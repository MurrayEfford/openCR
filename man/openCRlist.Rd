\name{openCRlist}
\alias{openCRlist}
\alias{[.openCRlist}

\title{Bundle openCR Models}

\description{
  
  Fitted models are bundled together for convenience.
  
}

\usage{

openCRlist (...)
\method{[}{openCRlist}(x, i)

}

\arguments{
  
  \item{\dots}{ openCR objects}
  \item{x}{openCRlist}
  \item{i}{indices}

}

\details{

\code{openCRlist} forms a special list (class `openCRlist') of fitted model (openCR) objects.
This may be used as an argument of \code{AIC}, \code{predict}, \code{make.table} etc.

Methods are provided for the generic function \code{c} and list extraction `['.

}

\value{

openCRlist object

}

\seealso{
  \code{\link{AIC.openCR}}
  \code{\link{predict.openCR}}
  \code{\link{make.table}}
}

\examples{

\dontrun{
fit0 <- openCR.fit (dipperCH)
fitt <- openCR.fit (dipperCH, model=phi~t)
fits <- openCRlist(fit0,fitt)
AIC(fits)
make.table(fits, 'phi')
}

}
\keyword{ model }

