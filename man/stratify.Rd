\name{stratify}
\alias{stratify}

\title{Stratify Capture-Recapture Data}

\description{

Arrange existing capthist data in stratified form.

}

\usage{
stratify(..., intervals = NULL, MoreArgs = list(), covariate = NULL, bytraps = FALSE)
}

\arguments{
  \item{\dots}{one or more multi-session capthist objects, or a list of such objects}
  \item{intervals}{list of intervals vectors, one for each multi-session capthist in \dots}
  \item{MoreArgs}{list of other arguments passed to \code{\link[secr]{join}}}
  \item{covariate}{character; name of individual or trap covariate to stratify by}
  \item{bytraps}{logical; if TRUE then covariate is interpreted as the name of a detector covariate}
}

\details{

 The argument \dots may be
\enumerate{
\item{a list of single-session capthist, one for each stratum (sessions already joined)}
\item{a list of multi-session capthist, one for each stratum (sessions will be joined)}
\item{one single-session capthist, to split by \code{covariate} (sessions already joined)}
\item{one multi-session capthist, to be joined as one then split by \code{covariate}}
}

Cases 1 and 2 result in one stratum for each component of the input list. 
Cases 3 and 4 result in one stratum for each level of \code{covariate}.

The result in Case 1 is identical to \code{MS.capthist(\dots)}.
 
The argument \code{intervals} refers to the intervals between primary sessions 
before joining (Cases 2,4 only) (see Examples).

MoreArgs may include the arguments remove.dupl.sites, tol, sites.by.name or drop.sites of \code{\link[secr]{join}}; these otherwise take their default values.

}

\value{

Multi-stratum (multi-session) capthist object for which each component has been `join'ed.

}

\seealso{

\code{\link[secr]{join}}, 
\code{\link[secr]{MS.capthist}}, 
\code{\link{openCR.fit}}, 
\code{\link{strata}}

}

\examples{

# FebpossumCH comprises 9 annual February sessions.
# The individual covariate 'sex' takes values 'F' and 'M', resulting in two strata.
# 'intervals' can be omitted as the default does the same job.
ch <- stratify(FebpossumCH, covariate = 'sex', intervals = rep(list(rep(1,8)),2))
summary(ch, terse = TRUE)

}

\keyword{ manip }
