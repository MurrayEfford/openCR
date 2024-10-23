\name{AIC.openCR}
\alias{AIC.openCR}
\alias{AIC.openCRlist}
\alias{logLik.openCR}

\title{ Compare openCR Models }

\description{

Terse report on the fit of one or more spatially explicit
capture--recapture models. Models with smaller values of AIC (Akaike's
Information Criterion) are preferred.

}

\usage{

\method{AIC}{openCR}(object, ..., sort = TRUE, k = 2, dmax = 10,  use.rank = FALSE,
                        svtol = 1e-5, criterion = c('AIC','AICc'), n = NULL)

\method{AIC}{openCRlist}(object, ..., sort = TRUE, k = 2, dmax = 10,  use.rank = FALSE,
                        svtol = 1e-5, criterion = c('AIC','AICc'), n = NULL)

\method{logLik}{openCR}(object, ...)

}

\arguments{
  \item{object}{ \code{openCR} object output from the function
      \code{\link{openCR.fit}}, or \code{openCRlist}}
  \item{\dots}{ other \code{openCR} objects }
  \item{sort}{ logical for whether rows should be sorted by ascending AICc }
  \item{k}{ numeric, the penalty per parameter to be used; always k = 2 in this method}
  \item{dmax}{ numeric, the maximum AIC difference for inclusion in
    confidence set}
  \item{use.rank}{ logical; if TRUE the number of parameters is based on
    the rank of the Hessian matrix}
  \item{svtol}{ minimum singular value (eigenvalue) of Hessian used when
    counting non-redundant parameters }
  \item{criterion}{ character, criterion to use for model comparison and
  weights}
  \item{n}{integer effective sample size}
}

\details{

Models to be compared must have been fitted to the same data and use the
same likelihood method (full vs conditional).

AIC with small sample adjustment is given by

\deqn{ \mbox{AIC}_c = -2\log(L(\hat{\theta})) + 2K +
\frac{2K(K+1)}{n-K-1} }{AICc = -2log(L(theta-hat)) + 2K +
2K(K+1)/(n-K-1)}

where \eqn{K} is the number of ``beta" parameters estimated. By default, the effective sample size \eqn{n} is the number of individuals observed at least once (i.e. the
number of rows in \code{capthist}). This differs from the default in MARK which for CJS models is the sum of the sizes of release cohorts (see \code{\link{m.array}}).

Model weights are calculated as \deqn{w_i = \frac{\exp(-\Delta_i / 2)}{
\sum{\exp(-\Delta_i / 2)}}}{w_i = exp(-dAICc_i / 2) / sum{ exp(-dAICc_i
/ 2) }}

Models for which dAIC > \code{dmax} are given a weight of zero and are
excluded from the summation. Model weights may be used to form
model-averaged estimates of real or beta parameters with
\code{\link{modelAverage}} (see also Buckland et al. 1997, Burnham and
Anderson 2002).

The argument \code{k} is included for consistency with the generic
method \code{AIC}.

}

\value{
  A data frame with one row per model.  By default, rows are sorted by ascending AIC.
  \item{model }{character string describing the fitted model}
  \item{npar }{number of parameters estimated}
  \item{rank }{rank of Hessian}
  \item{logLik }{maximized log likelihood}
  \item{AIC }{Akaike's Information Criterion}
  \item{AICc }{AIC with small-sample adjustment of Hurvich & Tsai (1989)}
  \item{dAICc }{difference between AICc of this model and the one with smallest AIC}
  \item{AICwt }{AICc model weight}

  \code{logLik.openCR} returns an object of class `logLik' that has
  attribute \code{df} (degrees of freedom = number of estimated
  parameters).
  
}
\references{ 

Buckland S. T., Burnham K. P. and Augustin, N. H. (1997) Model selection: an integral part of inference. \emph{Biometrics} \bold{53}, 603--618.

Burnham, K. P. and Anderson, D. R. (2002) \emph{Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach}. Second edition. New York: Springer-Verlag.

Hurvich, C. M. and Tsai, C. L. (1989) Regression and time series model selection in small samples. \emph{Biometrika} \bold{76}, 297--307.


}

\note{

The default criterion is AIC, not AICc as in \pkg{secr} 3.1.

Computed values differ from MARK for various reasons. MARK uses the
number of observations, not the number of capture histories when
computing AICc. It is also likely that MARK will count parameters
differently.
  
It is not be meaningful to compare models by AIC if they relate to
different data.

The issue of goodness-of-fit and possible adjustment of AIC for
overdispersion has yet to be addressed (cf QAIC in MARK).

}

\seealso{
  \code{\link[stats]{AIC}}, \code{\link{openCR.fit}},
  \code{\link{print.openCR}}, \code{\link[secr]{LR.test}}
}

\examples{

\dontrun{
m1 <- openCR.fit(ovenCH, type = 'JSSAf')
m2 <- openCR.fit(ovenCH, type = 'JSSAf', model = list(p~session))
AIC(m1, m2)
}

}

\keyword{ models }

