\name{print.derivedopenCR}
\alias{print.derivedopenCR}

\title{Print Method for Derived Estimates}

\description{

Formats output from \code{\link{derived.openCR}}.

}

\usage{

\method{print}{derivedopenCR}(x, Dscale = NULL, legend = FALSE, ...)

}

\arguments{
  \item{x}{object from \code{derived.openCR}}
  \item{Dscale}{numeric optional multiplier for densities (overrides saved Dscale)}
  \item{legend}{logical. if TRUE then a legend is provided to column headings}
  \item{\dots}{other arguments passed to \code{\link{print.data.frame}} }
}

\details{

By default (i.e. when not not specified in the in the \dots argument), 
\code{row.names = FALSE} and \code{digits = 4}.

}

\seealso{

\code{\link{derived.openCR}}

}

