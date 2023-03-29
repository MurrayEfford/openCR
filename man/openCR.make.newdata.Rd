\name{makeNewData}
\alias{makeNewData}
\alias{makeNewData.openCR}

\title{ Create Default Design Data }

\description{

Internal function used to generate a dataframe containing design data
for the base levels of all predictors in an openCR object.

}

\usage{

\method{makeNewData}{openCR}(object, all.levels = FALSE, ...)

}

\arguments{

\item{object}{fitted openCR model object}
\item{all.levels}{logical; if TRUE then all covariate factor levels appear in the output}
\item{\dots}{other arguments (not used)}

}

\details{

\code{makeNewData} is used by \code{predict} in lieu of
user-specified `newdata'. There is seldom any need to call
\code{makeNewData} directly.

\code{makeNewData} respects the details argument `agebreaks' for grouping ages 
(\pkg{openCR} >= 2.2.6).

}

\value{

A dataframe with one row for each session, and columns for the
predictors used by \code{object$model}.
  
}

\seealso{
  
\code{\link{openCR.fit}}
  
}

\examples{

\dontrun{

## null example (no covariates)
ovenCJS <- openCR.fit(ovenCH)
makeNewData(ovenCJS)
%openCR.make.newdata(ovenCJS)

}

}

\keyword{ models }
