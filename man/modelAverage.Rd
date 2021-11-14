\name{modelAverage}
\alias{modelAverage}
\alias{modelAverage.openCR}
\alias{modelAverage.openCRlist}

\title{ Averaging of OpenCR Models Using Akaike's Information Criterion}

\description{
  
AIC- or AICc-weighted average of estimated `real' or `beta' parameters
from multiple fitted openCR models. 

For a technical reason that will disappear when the next version of \pkg{secr} is 
released, these are internal functions are not exported by \pkg{openCR} and 
for now must be called with the prefix `openCR:::'.

}

\usage{

\method{modelAverage}{openCR}(object, ..., realnames = NULL, betanames = NULL, 
    newdata = NULL, alpha = 0.05, dmax = 10, covar = FALSE, average = c("link", 
    "real"), criterion = c("AICc","AIC"), CImethod = c("Wald", "MATA"))
    
\method{modelAverage}{openCRlist}(object, ..., realnames = NULL, betanames = NULL, 
    newdata = NULL, alpha = 0.05, dmax = 10, covar = FALSE, average = c("link", 
    "real"), criterion = c("AICc","AIC"), CImethod = c("Wald", "MATA"))

}

\arguments{
\item{object}{\code{openCR} or \code{openCRlist} objects}
  \item{\dots}{ other \code{openCR} objects (modelAverage.openCR() only)}
  \item{realnames}{ character vector of real parameter names }
  \item{betanames}{ character vector of beta parameter names }
  \item{newdata}{ optional dataframe of values at which to evaluate models }
  \item{alpha}{ alpha level for confidence intervals }
  \item{dmax}{ numeric, the maximum AIC or AICc difference for inclusion in confidence set }
  \item{covar}{ logical, if TRUE then return variance-covariance matrix }
  \item{average}{ character string for scale on which to average real
    parameters }
  \item{criterion}{ character, information criterion to use for model weights}  
  \item{CImethod}{ character, type of confidence interval (see Details)}
}

\details{

Models to be compared must have been fitted to the same data and use the
same likelihood method (full vs conditional). If \code{realnames} =
NULL and \code{betanames} = NULL then all real parameters will be
averaged; in this case all models must use the same real parameters. To
average beta parameters, specify \code{betanames} (this is ignored if a
value is provided for \code{realnames}). See \code{\link{predict.openCR}}
for an explanation of the optional argument \code{newdata};
\code{newdata} is ignored when averaging beta parameters. 

Model-averaged estimates for parameter \eqn{\theta}{theta} are given by
\deqn{\hat{\theta} = \sum\limits _k w_k \hat{\theta}_k}{theta-hat = sum(
w_k * theta-hat_k)} where the subscript \eqn{k} refers to a specific
model and the \eqn{w_k} are AIC or AICc weights (see
\code{\link{AIC.openCR}} for details). Averaging of real parameters may be
done on the link scale before back-transformation
(\code{average="link"}) or after back-transformation
(\code{average="real"}).

Models for which dAIC > \code{dmax} (or dAICc > \code{dmax}) are given a
weight of zero and effectively are excluded from averaging.

Also, \deqn{\mbox{var} (\hat{\theta}) = \sum\limits _{k} { w_{k}
( \mbox{var}(\hat{\theta}_{k} | \beta _k) + \beta _k ^2)}
}{var(theta-hat) = sum(w_k (var(theta-hat_k) + beta_k^2))}

where \eqn{\hat{\beta} _k = \hat{\theta}_k - \hat{\theta}}{beta-hat_k =
theta-hat_k -- theta-hat} and the variances are asymptotic estimates
from fitting each model \eqn{k}. This follows Burnham and Anderson
(2004) rather than Buckland et al. (1997).

Two methods are offered for confidence intervals. The default `Wald'
uses the above estimate of variance. The alternative `MATA'
(model-averaged tail area) avoids estimating a weighted variance and
is thought to provide better coverage at little cost in increased
interval length (Turek and Fletcher 2012). Turek and Fletcher (2012)
also found averaging with AIC weights (here \code{criterion = 'AIC'})
preferable to using AICc weights, even for small
samples. \code{CImethod} does not affect the reported standard errors.

}

\value{

A list (one component per parameter) of model-averaged estimates, their
standard errors, and a \eqn{100(1-\alpha)}{100(1-alpha)}\% confidence
interval. The interval for real parameters is backtransformed from the
link scale. If there is only one row in \code{newdata} or beta
parameters are averaged or averaging is requested for only one parameter
then the array is collapsed to a matrix. If \code{covar = TRUE} then a
list is returned with separate components for the estimates and the
variance-covariance matrices.

}

\references{ 

Buckland S. T., Burnham K. P. and Augustin, N. H. (1997) Model
selection: an integral part of inference. \emph{Biometrics} \bold{53},
603--618.

Burnham, K. P. and Anderson, D. R. (2002) \emph{Model Selection and
Multimodel Inference: A Practical Information-Theoretic Approach}.
Second edition. New York: Springer-Verlag.

Burnham, K. P. and Anderson, D. R. (2004) Multimodel inference -
understanding AIC and BIC in model selection. \emph{Sociological Methods
& Research} \bold{33}, 261--304.

Turek, D. and Fletcher, D. (2012) Model-averaged Wald confidence
intervals. \emph{Computational statistics and data analysis} \bold{56},
2809--2815.

}

\seealso{ 

\code{\link{AIC.openCR}}, 
\code{\link{make.table}}, 
\code{\link{openCR.fit}}, 
\code{\link{openCRlist}}

}

\examples{

## Compare two models fitted previously

cjs1 <- openCR.fit(dipperCH, model=p~1)
cjs2 <- openCR.fit(dipperCH, model=p~session)
openCR:::modelAverage(cjs1, cjs2)

## or
cjs12 <- openCRlist(cjs1, cjs2)
openCR:::modelAverage(cjs12)

}

\keyword{ models }

