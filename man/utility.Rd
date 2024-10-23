\name{JS.counts}
\alias{JS.counts}
\alias{m.array}
\alias{bd.array}

\title{Summarise Non-spatial Open-population Data}

\description{

Simple conventional summaries of data held in secr `capthist' objects.

}

\usage{

JS.counts(object, primary.only = TRUE, stratified = FALSE)
m.array(object, primary.only = TRUE, never.recaptured = TRUE, 
    last.session = TRUE, stratified = FALSE)
bd.array(beta, phi)

}

\arguments{
  \item{object}{secr capthist object or similar}
  \item{primary.only}{logical; if TRUE then counts are tabuated for primary sessions }
  \item{stratified}{logical; if TRUE then sessions of multisession object summarised separately }
  \item{never.recaptured}{logical; if TRUE then a column is added for animals never recaptured }
  \item{last.session}{logical; if TRUE releases are reported for the last session}
  \item{beta}{numeric vector of entry probabilities, one per primary session}
  \item{phi}{numeric vector of survival probabilities, one per primary session}

}

\details{

The input is a capthist object representing a multi-session capture--recapture study. This may be (i) a single-session capthist in which occasions are understood to represent primary sessions, or (ii) a multi-session capthist object that is automatically converted to a single session object with \code{\link[secr]{join}} (any secondary sessions (occasions) are first collapsed with \code{reduce(object, by = 'all')}*, or (iii) a multi-session capthist object in which sessions are interpreted as strata.

The argument \code{primary.only} applies for single-session input with a robust-design structure defined by the \code{\link[secr]{intervals}}. \code{last.session} results in a final row with no recaptures.

If the covariates attribute of \code{object} includes a column named `freq' then this is used to expand the capture histories.

Conventional Jolly--Seber estimates may be computed with \code{\link{JS.direct}}.

\code{bd.array} computes the probability of each possible combination of birth and death times (strictly, the primary session at which an animal was first and last available for detection), given the parameter vectors beta and phi. These cell probabilities are integral to JSSA models.

* this may fail with nonspatial data.

}

\value{
For \code{JS.counts}, a data.frame where rows correspond to sessions and columns hold counts as follows --
\item{n}{number of individuals detected}
\item{R}{number of individuals released}
\item{m}{number of previously marked individuals}
\item{r}{number of released individuals detected in later sessions}
\item{z}{number known to be alive (detected before and after) but not detected in current session}

For \code{m.array}, a table object with rows corresponding to release cohorts and columns  corresponding to first--recapture sessions. The size of the release cohort is shown in the first column. Cells in the lower triangle have value NA and print as blank by default.

}

\seealso{

\code{\link[secr]{join}}, \code{\link{JS.direct}}

}

\examples{

JS.counts(ovenCH)
m.array(ovenCH)

\dontrun{

## probabilities of b,d pairs
fit <- openCR.fit(ovenCH, type = 'JSSAbCL')
beta <- predict(fit)$b$estimate
phi <- predict(fit)$phi$estimate
bd.array(beta, phi)

}

}
\keyword{ manip }