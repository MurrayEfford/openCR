\name{predict.openCR}
\alias{predict.openCR}
\alias{predict.openCRlist}

\title{ openCR Model Predictions }

\description{

Evaluate an openCR capture--recapture model. That is, compute the `real' parameters corresponding to the `beta' parameters of a fitted model for arbitrary levels of any variables in the linear predictor.

}

\usage{

\method{predict}{openCR} (object, newdata = NULL, se.fit = TRUE, alpha = 0.05, savenew = FALSE, ...) 
\method{predict}{openCRlist} (object, newdata = NULL, se.fit = TRUE, alpha = 0.05, savenew = FALSE, ...) 

}

\arguments{
  \item{object}{ \code{openCR} object output from \code{openCR.fit} }
  \item{newdata}{ optional dataframe of values at which to evaluate model }
  \item{se.fit}{ logical for whether output should include SE and confidence intervals }
  \item{alpha}{ alpha level }
  \item{savenew}{ logical; if TRUE then newdata is saved as an attribute }
  \item{\dots}{ other arguments passed to \code{\link{makeNewData}} }
  %\item{\dots}{ other arguments passed to \code{\link{openCR.make.newdata}} }
 }

\details{

Predictions are provided for each row in `newdata'. The default (constructed by 
\code{\link{makeNewData}}) is to limit those rows to the first-used level of 
factor predictors; to include all levels pass \code{all.levels = TRUE} to 
\code{makeNewData} in the \dots argument.

%Predictions are provided for each row in `newdata'. The default (constructed by 
%\code{\link{openCR.make.newdata}}) is to limit those rows to the first-used level of 
%factor predictors; to include all levels pass \code{all.levels = TRUE} to 
%\code{openCR.make.newdata} in the \dots argument.

}

\seealso{ \code{\link{AIC.openCR}}, \code{\link{openCR.fit}} }

\examples{

\dontrun{

c1 <- openCR.fit(ovenCH, type='CJS', model=phi~session)
predict(c1)

}

}

