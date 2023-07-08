\name{make.table}
\alias{make.table}

\title{Tabulate Estimates From Multiple Models}

\description{

Session-specific estimates of real parameters (p, phi, etc.) are arranged in a rectangular table.

}

\usage{

make.table(fits, parm = "phi", fields = "estimate", strata = 1, 
    collapse = FALSE, ...)

}

\arguments{
  \item{fits}{openCRlist object}
  \item{parm}{character name of real parameter estimate to tabulate}
  \item{fields}{character column from predict (estimate, SE.estimate, lcl, ucl)}
  \item{strata}{integer; indices of strata to report}
  \item{collapse}{logical; if TRUE stratum-specific results are collapsed to single table}
  \item{...}{arguments passed to \code{\link{predict.openCRlist}}}
}

\details{

The input will usually be from \code{par.openCR.fit}.

\code{\link{collate.openCR}} is a flexible alternative.

}

\value{

A table object.

}

\seealso{

\code{\link{collate.openCR}},
\code{\link{par.openCR.fit}}, 
\code{\link{openCRlist}}

}

\examples{

\dontrun{

arglist <- list(
    constant = list(capthist = ovenCHp, model = phi~1), 
    session.specific = list(capthist = ovenCHp, model = phi~session)
)
fits <- par.openCR.fit(arglist, trace = FALSE)
print(make.table(fits), na = ".")

}

}

\keyword{ manip }