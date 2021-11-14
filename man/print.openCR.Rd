\name{print.openCR}
\alias{print.openCR}
\alias{summary.openCR}

\title{ Print or Summarise openCR Object }

\description{

  Print results from fitting a spatially explicit capture--recapture model, or generate a list of summary data.

}

\usage{
\method{print}{openCR} (x, newdata = NULL, alpha = 0.05, svtol = 1e-5,...)
\method{summary}{openCR} (object, newdata = NULL, alpha = 0.05, svtol = 1e-5, deriv = FALSE, ...)
}

\arguments{
  \item{x}{ \code{openCR} object output from \code{openCR.fit} }
  \item{object}{ \code{openCR} object output from \code{openCR.fit} }
  \item{newdata}{ optional dataframe of values at which to evaluate model }
  \item{alpha}{ alpha level }
  \item{svtol}{ threshold for non-null eigenvalues when computing numerical rank }
  \item{deriv}{logical; if TRUE then table of derived parameters is calculated}
  \item{...}{ other arguments passed to \code{\link{derived.openCR}} by \code{summary.openCR}}
}

\details{

Results are potentially complex and depend upon the analysis (see below). Optional \code{newdata} should be a dataframe with a column for each of the variables in the model. If \code{newdata} is missing then a dataframe is constructed automatically.  Default \code{newdata} are for a naive animal on the first occasion; numeric covariates are set to zero and factor covariates to their base (first) level. Confidence intervals are 100 (1 -- alpha) \% intervals.

\tabular{ll}{
call           \tab the function call \cr
time           \tab date and time fitting started \cr
N animals      \tab number of distinct animals detected \cr
N captures     \tab number of detections \cr
N sessions     \tab number of sampling occasions \cr
Model          \tab model formula for each `real' parameter \cr
Fixed          \tab fixed real parameters \cr
N parameters   \tab number of parameters estimated  \cr
Log likelihood \tab log likelihood \cr
AIC            \tab Akaike's information criterion \cr
AICc           \tab AIC with small sample adjustment (Burnham and Anderson 2002) \cr
Beta parameters \tab coef of the fitted model, SE and confidence
intervals \cr
Eigenvalues     \tab scaled eigenvalues of Hessian matrix (maximum 1.0) \cr
Numerical rank  \tab number of eigenvalues exceeding svtol \cr
vcov            \tab variance-covariance matrix of beta parameters  \cr
Real parameters \tab fitted (real) parameters evaluated at base levels of covariates \cr
}

AICc is computed with the default sample size (number of individuals) and parameter count (use.rank = FALSE).

}


\value{

The \code{summary} method constructs a list of outputs similar to those printed by the \code{print} method, 
but somewhat more concise and re-usable:

\tabular{ll}{
versiontime     \tab secr version, and date and time fitting started \cr
traps*           \tab detector summary \cr
capthist        \tab capthist summary (primary and secondary sessions, numbers of animals and detections)\cr
intervals        \tab intervals between primary sessions \cr
mask*            \tab mask summary \cr
modeldetails    \tab miscellaneous model characteristics (type etc.) \cr
AICtable        \tab single-line output of AIC.openCR \cr
coef            \tab table of fitted coefficients with CI \cr
predicted       \tab predicted values (`real' parameter estimates) \cr
derived         \tab output of derived.openCR (optional) \cr
}
* spatial models only
}

\references{ 
Burnham, K. P. and Anderson, D. R. (2002) \emph{Model selection and multimodel inference: a practical information-theoretic approach}. Second edition. New York: Springer-Verlag.
}

\seealso{ \code{\link{AIC.openCR}}, \code{\link{openCR.fit}} }

\examples{

\dontrun{

c1 <- openCR.fit(ovenCH, type='CJS', model=phi~session)
c1

}

}

\keyword{ print }
