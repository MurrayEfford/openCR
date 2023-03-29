\name{openCR.design}
\alias{openCR.design}

\title{Design Data for Open population Models}

\description{
  Internal function used by \code{\link{openCR.fit}}. 
}
\usage{

openCR.design(capthist, models, type, naive = FALSE, stratumcov = NULL, 
    sessioncov = NULL, timecov = NULL, agecov = NULL, 
    dframe = NULL, contrasts = NULL, initialage = 0, 
    minimumage = 0, maximumage = 1, agebreaks = NULL, CJSp1 = FALSE, ...)

}

\arguments{
  \item{capthist}{ single-session \code{capthist} object }
  \item{models}{ list of formulae for parameters of detection }
  \item{type}{character string for type of analysis "CJS", "JSSAfCL" etc. (see \code{\link{openCR.fit}})}
  \item{naive}{ logical if TRUE then modelled parameter is for a naive
    animal (not caught previously)}
  \item{timecov}{ optional vector or dataframe of values of occasion-specific covariate(s). }
  \item{stratumcov}{ optional dataframe of values of stratum-specific covariate(s) }
  \item{sessioncov}{ optional dataframe of values of session-specific covariate(s) }
  \item{agecov}{ optional dataframe of values of age-specific covariate(s) }
  
  \item{dframe}{ optional data frame of design data for detection parameters }
  \item{contrasts}{contrast specification as for \code{\link[stats]{model.matrix}}}
  \item{initialage}{numeric or character (name of individual covariate containing initial ages)}
  \item{minimumage}{numeric; ages younger than minimum are truncated up}
  \item{maximumage}{numeric; ages older than maximum are truncated down}
  \item{agebreaks}{numeric vector of age-class limits}
  \item{CJSp1}{logical; if TRUE detection is modelled on first primary session in CJS models}
  \item{\dots}{ other arguments passed to the \R function \code{\link{model.matrix}} }
}

\details{

  This is an internal \pkg{openCR} function that you are unlikely ever
to use. \dots may be used to pass \code{contrasts.arg} to
\code{model.matrix}.

Each real parameter is notionally different for each unique combination
of individual, secondary session, detector and latent class, i.e., for \eqn{n} individuals,
\eqn{S} secondary sessions, \eqn{K} detectors and \eqn{m} latent classes there are
\emph{potentially} \eqn{n \times S \times K \times m }{n x S x K x m} different
values. Actual models always predict a much reduced set of distinct
values, and the number of rows in the design matrix is reduced
correspondingly; a parameter index array allows these to retrieved for
any combination of individual, session and detector. 

\code{openCR.design} is less tolerant than \code{openCR.fit} regarding
the inputs `capthist' and `models'. Model formulae are processed by \code{openCR.fit}
to a standard form (a named list of formulae) before they are passed to
\code{openCR.design}, and multi-session capthist objects are
automatically `reduced' and `joined' for open-population analysis.

If \code{timecov} is a single vector of values (one for each secondary session) 
then it is treated as a covariate named `tcov'.
If \code{sessioncov} is a single vector of values (one for each primary session) 
then it is treated as a covariate named `scov'.

The \code{initialage} and \code{maximumage} arguments are usually passed via the
\code{openCR.fit} `details' argument.

\code{agecov} may be used to group ages. It should have length (or number of rows) 
equal to \code{maximumage} + 1. Alternatively, age classes may be defined with the argument \code{agebreaks}; this is preferred from openCR 2.2.6.

}

\value{
  
  A list with the components
  
  \item{designMatrices}{ list of reduced design matrices, one for each
    real parameter }
  
  \item{parameterTable}{ index to row of the reduced design matrix for
  each real parameter; dim(parameterTable) = c(uniquepar, np),
  where uniquepar is the number of unique combinations of paramater
  values (uniquepar < \eqn{nSKM}) and np is the number of parameters in
  the detection model. }

  \item{PIA}{ Parameter Index Array - index to row of parameterTable for
    a given animal, occasion and latent class; dim(PIA) = c(n,S,K,M) }

  \item{validlevels}{ for J primary sessions, a logical matrix of np rows and J columns, mostly
    TRUE, but FALSE for impossible combinations e.g. CJS recapture
    probability in session 1 (validlevels["p",1]) unless \code{CJSp1 = TRUE}, or CJS final survival
    probability (validlevels["phi",J]). Also, validlevels["b",1] is
    FALSE with type = "JSSA..." because of the constraint that entry
    parameters sum to one. }
    
  \item{individual}{ TRUE if uses individual variate(s)}  
  
  \item{agelevels}{levels for age factor (cut numeric ages) if `age' in model}

  }

  \note{The component validlevels is TRUE in many cases for which a
    parameter is redundant or confounded (e.g. validlevels["phi",J-1]);
    these are sorted out `post hoc' by examining the fitted values,
    their asymptotic variances and the eigenvalues of the Hessian
    matrix.}

  \seealso{

  \code{\link{openCR.fit}}
  
}

\examples{

## this happens automatically in openCR.fit
ovenCH1 <- join(reduce(ovenCH, by = "all", newtraps=list(1:44)))

openCR.design (ovenCH1, models = list(p = ~1, phi = ~session),
    interval = c(1,1,1,1), type = "CJS")

}

\keyword{ manip }
