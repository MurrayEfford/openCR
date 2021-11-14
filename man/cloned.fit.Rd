\name{cloned.fit}
\alias{cloned.fit}

\title{

Cloning to Evaluate Identifiability

}

\description{

The identifiability of parameters may be examined by refitting a model with cloned data (each capture history replicated \code{nclone} times). For identifiable parameters the estimated variances are proportional to \code{1/nclone}.

}

\usage{

cloned.fit(object, nclone = 100, newdata = NULL, linkscale = FALSE)

}

\arguments{
  \item{object}{previously fitted openCR object}
  \item{nclone}{integer number of times to replicate each capture history}
  \item{newdata}{optional dataframe of values at which to evaluate model }
  \item{linkscale}{logical; if TRUE then comparison uses SE of linear predictors}
}

\details{

The key output is the ratio of SE for estimates from the uncloned and cloned datasets, adjusted for the level of cloning (\code{nclone}). For identifiable parameters the ratio is expected to be 1.0.

Cloning is not implemented for spatial models.

The comparison may be done either on the untransformed scale (using approximate SE) or on the link scale.

}
\value{

Dataframe with columns* --
  \item{estimate }{original estimate }
  \item{SE.estimate }{original SE}
  \item{estimate.xxx}{cloned estimate (xxx = nclone)}
  \item{SE.estimate.xxx }{cloned SE}
  \item{SE.ratio}{SE.estimate / SE.estimate.xxx / sqrt(nclone)}
  
* `estimate' becomes `beta' when \code{linkscale = TRUE}.

}

\references{

Lele, S.R., Nadeem, K. and Schmuland, B. (2010) Estimability and likelihood inference for generalized
linear mixed models using data cloning. \emph{Journal of the American Statistical Association} \bold{105}, 1617--1625.

}

\seealso{

\code{\link{openCR.fit}}

}

\examples{

\dontrun{
fit <- openCR.fit(dipperCH)
cloned.fit(fit)
}

}

\keyword{ model }