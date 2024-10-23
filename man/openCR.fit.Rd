\name{openCR.fit}
\alias{openCR.fit}

\title{ Fit Open Population Capture--Recapture Model }

\description{

  Nonspatial or spatial open-population analyses are performed on data
formatted for `secr'. Several parameterisations are provided for the
nonspatial Jolly-Seber Schwarz-Arnason model (`JSSA', also known as
`POPAN'). Corresponding spatial models are designated
`JSSAsecr'. The prefix `PLB' (Pradel-Link-Barker) is used for versions of the 
JSSA models that are conditional on the number observed. Cormack-Jolly-Seber 
(CJS) models are also fitted.

}

\usage{

openCR.fit (capthist, type = "CJS", model = list(p~1, phi~1, sigma~1), 
    distribution = c("poisson", "binomial"), mask = NULL, 
    detectfn = c("HHN", "HHR", "HEX", "HAN", "HCG", "HVP", "HPX"), binomN = 0, 
    movementmodel = c('static', 'BVN', 'BVE', 'BVT', 'RDE', 'RDG','RDL','IND', 'UNI',
      'BVNzi', 'BVEzi', 'RDEzi', 'INDzi', 'UNIzi'), edgemethod = 
    c("truncate", "wrap", "none"), kernelradius = 30, sparsekernel = TRUE, 
    start = NULL, link = list(), fixed = list(), stratumcov = NULL, 
    sessioncov = NULL, timecov = NULL, agecov = NULL, dframe = NULL, 
    dframe0 = NULL, details = list(), method = "Newton-Raphson", trace = NULL, 
    ncores = NULL, stratified = FALSE, ...)
}
% changed 2.1.1 2021-10-19
% c("truncate","wrap","none"), kernelradius = 10, sparsekernel = FALSE, 

\arguments{

  \item{capthist}{ \code{capthist} object from `secr'}

  \item{type}{character string for type of analysis (see Details)}

  \item{model}{ list with optional components, each symbolically
  defining a linear predictor for the relevant real parameter using
  \code{formula} notation. See Details for names of real parameters. }
  
  \item{distribution}{character distribution of number of individuals detected}

  \item{mask}{ single-session \code{\link[secr]{mask}} object; required for spatial (secr) models }

  \item{detectfn}{character code}

  \item{binomN}{ integer code for distribution of counts (see \code{\link[secr]{secr.fit}}) }

  \item{movementmodel}{character; model for movement between primary sessions (see Details) }
  
  \item{edgemethod}{character; method for movement at edge of mask (see Details)}
  
  \item{kernelradius}{integer; radius in mask cells of discretized kernel (movement models only)}
  
  \item{sparsekernel}{logical; if TRUE then only cardinal and intercardinal axes are included}

  \item{start}{ vector of initial values for beta parameters, or fitted
    model(s) from which they may be derived }
  
  \item{link}{ list with named components, each a character string in
  \{"log", "logit", "loglog", "identity", "sin", "mlogit"\} for the link function
  of the relevant real parameter }

 \item{fixed}{ list with optional components corresponding to each
 `real' parameter, the scalar value to which parameter is to be fixed }

  \item{stratumcov}{ optional dataframe of values of stratum-specific
  covariate(s). }
  
  \item{sessioncov}{ optional dataframe of values of session-specific
  covariate(s). }
  
  \item{timecov}{ optional dataframe of values of occasion-specific
  covariate(s). }

  \item{agecov}{ optional dataframe of values of age-specific covariate(s) }

\item{dframe}{ optional data frame of design data for detection
  parameters (seldom used) }

\item{dframe0}{ optional data frame of design data for detection
  parameters of naive (undetected) animals (seldom used) }

  \item{details}{ list of additional settings (see Details) }

  \item{method}{ character string giving method for maximizing log
    likelihood }
  
  \item{trace}{ logical or integer; output log likelihood at each evaluation, 
  or at some lesser frequency as given}

  \item{ncores}{integer number of cores for parallel processing (see Details) }
  
  \item{stratified}{logical; if TRUE then sessions of capthist interpreted 
  as indpendent strata}
  
  \item{\dots}{ other arguments passed to join() }
}

\details{

The permitted nonspatial models are CJS, Pradel, Pradelg, JSSAbCL = PLBb, JSSAfCL = PLBf, JSSAgCL = PLBg, JSSAlCL = PLBl, JSSAb, JSSAf, JSSAg, JSSAl, JSSAB and JSSAN. 

The permitted spatial models are CJSsecr, JSSAsecrbCL = PLBsecrb, JSSAsecrfCL = PLBsecrf, JSSAsecrgCL = PLBsecrg, JSSAsecrlCL = PLBsecrl, JSSAsecrb, JSSAsecrf, JSSAsecrg, JSSAsecrl, JSSAsecrB, JSSAsecrN, secrCL, and secrD. 

See \href{../doc/openCR-vignette.pdf}{openCR-vignette.pdf} for a table of the `real' parameters associated with each model type.

Parameterisations of the JSSA models differ in how they include
recruitment: the core parameterisations express recruitment either as a
per capita rate (`f'), as a finite rate of increase for the population
(`l' for lambda) or as per-occasion entry probability (`b' for the
classic JSSA beta parameter, aka PENT in MARK). Each of these models may
be fitted by maximising either the full likelihood, or the likelihood
conditional on capture in the Huggins (1989) sense, distinguished by the
suffix `CL'. Full-likelihood JSSA models may also be parameterized in
terms of the time-specific absolute recruitment (BN, BD) or the
time-specific population size(N) or density (D).

`secrCL' and `secrD' are closed-population spatial models.

Data are provided as \pkg{secr} `capthist' objects, with some
restrictions. For nonspatial analyses, `capthist' may be
single-session or multi-session, with any of the main detector types. For
spatial analyses `capthist' should be a single-session dataset of a point 
\link[secr]{detector} type (`multi', `proximity' or `count') (see also
details$distribution below). In openCR the occasions of a single-session
dataset are treated as open-population temporal samples except that occasions 
separated by an interval of zero (0) are from the same primary session (multi-session
input is collapsed to single-session if necessary).

\code{model} formulae may include the pre-defined terms
`session',`Session', `h2', and `h3' as in \pkg{secr}. `session'
is the name given to primary sampling times in `secr', so a fully
time-specific CJS model is \code{list(p ~ session, phi
~ session)}. `t' is a synonym of `session'. `Session' is for a
trend over sessions. `h2' and `h3' allow finite mixture models. 

Learned (behavioural) responses (`b', `B', etc.) were redefined and extended in 
version 1.3.0. The \href{../doc/openCR-vignette.pdf}{vignette} should be consulted 
for current definitions.

Formulae may also include named occasion-specific and session-specific 
covariates in the dataframe arguments `timecov' and `sessioncov' 
(occasion = secondary session of robust design). Named age-specific covariates 
in 'agecov` are treated similarly. Individual covariates present as an attribute of
the `capthist' input may be used in CJS and ..CL models. Groups are not
supported in this version, but may be implemented via a factor-level
covariate in ..CL models.

\code{distribution} specifies the distribution of the number of
individuals detected; this may be conditional on the population size (or number in the
masked area) ("binomial") or unconditional ("poisson").
\code{distribution} affects the sampling variance of the estimated
density. The default is "poisson" as in \pkg{secr}.

Movement models are list at \link{Movement models}. Their use is described in the \href{../doc/openCR-vignette.pdf}{vignette}. 

\code{edgemethod} controls movement probabilities at the mask edge in spatial models that include movement. "none" typically causes bias in estimates; "wrap" wraps kernel probabilities to the opposing edge of a rectangular mask; "truncate" scales the values of an edge-truncated kernel so that they always sum to 1.0 (safer and more general than "wrap").

The mlogit link function is used for the JSSA (POPAN) entry parameter 
`b' (PENT in MARK) and for mixture proportions, regardless of \code{link}.

Spatial models use one of the hazard-based detection functions (see \code{\link[secr]{detectfn}}) and require data
from independent point detectors (\pkg{secr} detector types `multi', `proximity' or `count').

Code is executed in multiple threads unless the user specifies \code{ncores = 1} or there is only one core available or \code{details$R == TRUE}. Setting \code{ncores = NULL} uses the existing value from the environment variable RCPP_PARALLEL_NUM_THREADS (see \code{\link[secr]{setNumThreads}}) or 2 if that has not been set. 

Optional stratification was introduced in \pkg{openCR} 2.0.0. See \href{../doc/openCR-vignette.pdf}{openCR-vignette.pdf} for details.

The \dots argument may be used to pass a vector of unequal intervals to 
join (\code{interval}), or to vary the tolerance for merging detector sites (\code{tol}).

The \code{start} argument may be 
\describe{
\item{- a vector of beta parameter values, one for each of the NP beta parameters in the model}{}
\item{- a named vector of beta parameter values in any order}{}
\item{- a named list of one or more real parameter values}{}
\item{- a single fitted secr or openCR model whose real parameters overlap with the current model}{}
\item{- a list of two fitted models}{}
}

In the case of two fitted models, the values are melded. This is handy for initialising an 
open spatial model from a closed spatial model and an open non-spatial model. If a beta 
parameter appears in both models then the first is used.

\code{details} is a list used for various specialized settings --

\tabular{lll}{
Component \tab Default \tab Description \cr

\code{agebreaks} \tab minimumage:maximumage \tab Limits of age classes (vector passed to \code{\link{cut}}) \cr

\code{autoini} \tab 1 \tab Number of the session used to determine 
initial values of D, lambda0 and sigma (secr types only) \cr

\code{CJSp1} \tab FALSE \tab Modified CJS model including initial detection 
(estimable with robust design and many spatial models) \cr

\code{contrasts} \tab NULL \tab Value suitable for the `contrasts.arg' argument 
of \code{\link{model.matrix}} used to specify the coding of factor predictors \cr

\code{control} \tab list() \tab Components may be named arguments of \code{nlm}, or passed intact as argument `control' of \code{optim} - useful for increasing maxit for \code{method = Nelder-Mead} (see vignette) \cr

\code{debug} \tab 0 \tab debug=1 prints various intermediate values; debug>=2 interrupts execution by calling browser() (position variable) \cr

\code{fixedbeta} \tab NULL \tab Vector with one element for each coefficient (beta parameter) in the model. Only 'NA' coefficients will be estimated; others will be fixed at the value given (coefficients define a linear predictor on the link scale). The number and order of coefficients may be determined by calling \code{openCR.fit} with trace = TRUE and interrupting execution after the first likelihood evaluation. \cr

\code{grain} \tab 1 \tab Obscure setting for multithreading - see \pkg{RcppParallel} package \cr
\code{hessian} \tab "auto" \tab Computation of the Hessian matrix from which 
variances and covariances are obtained. Options are "none" (no variances), 
"auto" or "fdhess" (use the function fdHess in \pkg{nlme}).  If "auto" then 
the Hessian from the optimisation function is used. \cr

\code{ignoreusage} \tab FALSE \tab Overrides usage in traps object of capthist \cr

\code{initialage} \tab 0 \tab Numeric (uniform age at first capture) 
or character value naming an individual covariate; see \code{\link{age.matrix}} \cr

\code{initialstratum} \tab 1 \tab Number of stratum to use for finding default starting values (cf autoini in \pkg{secr}) \cr

\code{LLonly} \tab FALSE \tab TRUE causes the function to return a single
evaluation of the log likelihood at the initial values, followed by the 
initial values \cr

\code{minimumage} \tab 0 \tab Sets a minimum age; see \code{\link{age.matrix}} \cr

\code{maximumage} \tab 1 \tab Sets a maximum age; older animals are recycled into 
this age class; see \code{\link{age.matrix}} \cr

\code{multinom} \tab FALSE \tab Include the multinomial constant in the 
reported log-likelihood. \cr

\code{r0} \tab 0.5 \tab effective radius of zero cell in movement kernel (multiple of cell width) \cr

\code{R} \tab FALSE \tab Switch from the default C++ code to 
slower functions in native R (useful for debugging; not all models) \cr

\code{squeeze} \tab TRUE \tab  Apply \code{\link{squeeze}} to capthist before analysis. Non-spatial models fit faster, because histories often non-unique.\cr

\code{userdist} \tab NULL \tab Function to compute distances (see \pkg{secr}) \cr

\code{stepmax} \tab NULL \tab stepmax argument of \code{\link[stats]{nlm}} (step on link scale) \cr

}

If \code{method = "Newton-Raphson"} then \code{\link[stats]{nlm}} is
used to maximize the log likelihood (minimize the negative log
likelihood); otherwise \code{\link[stats]{optim}} is used with the
chosen method ("BFGS", "Nelder-Mead", etc.).  If maximization fails a
warning is given appropriate to the method. \code{method = "none"} may 
be used to compute or re-compute the variance-covariance matrix at 
given starting values (i.e. providing a previously fitted model as 
the value of \code{start}).

Parameter redundancies are common in open-population models. The output
from \code{openCR.fit} includes the singular values (eigenvalues) of the
Hessian - a useful post-hoc indicator of redundancy (e.g., Gimenez et
al. 2004). Eigenvalues are scaled so the largest is 1.0. Very small
scaled values represent redundant parameters - in my experience with
simple JSSA models a threshold of 0.00001 seems effective.

[There is an undocumented option to fix specific `beta' parameters.]

Numeric ages may be grouped into age classes by providing `agebreaks'. 
In models, ~age then refers to the age-class factor. See the 
\href{../doc/openCR-vignette.pdf}{vignette} for more detail.

}

\value{
  
If \code{details$LLonly == TRUE} then a numeric vector is returned with logLik in 
position 1, followed by the named coefficients.

Otherwise, an object of class `openCR' with components

  \item{call }{function call }
  \item{capthist }{saved input (unique histories; see covariates(capthist)$freq for frequencies)}
  \item{type }{saved input}
  \item{model }{saved input}
  \item{distribution }{saved input}
  \item{mask }{saved input}
  \item{detectfn }{saved input}
  \item{binomN }{saved input}
  \item{movementmodel }{saved input}
  \item{edgemethod }{saved input}
  \item{usermodel }{saved input}
  \item{moveargsi }{relative positions of move.a and move.b arguments}
  \item{kernel}{coordinates of kernel (movement models only)}
  \item{start }{vector of starting values for beta parameters} 
  \item{link }{saved input}
  \item{fixed }{saved input}  
  \item{timecov }{saved input}
  \item{sessioncov }{saved input}
  \item{agecov }{saved input}
  \item{dframe }{saved input}
  \item{dframe0 }{saved input}
  \item{details }{saved input}
  \item{method }{saved input}
  \item{ncores }{saved input (NULL replaced with default)}
  \item{design }{reduced design matrices, parameter table and parameter
    index array for actual animals (see \code{\link{openCR.design}})}
  \item{design0 }{reduced design matrices, parameter table and parameter
    index array for `naive' animal (see \code{\link{openCR.design}})}
  \item{parindx }{list with one component for each real parameter giving
    the indices of the `beta' parameters associated with each real
    parameter}  
  \item{primaryintervals}{intervals between primary sessions}
  \item{vars }{vector of unique variable names in \code{model} }
  \item{betanames }{names of beta parameters}
  \item{realnames }{names of fitted (real) parameters }
  \item{sessionlabels}{name of each primary session}
  \item{fit}{list describing the fit (output from \code{nlm} or
    \code{optim}) }
  \item{beta.vcv }{variance-covariance matrix of beta parameters }  
  \item{eigH }{vector of eigenvalue corresponding to each beta parameter }
  \item{version }{openCR version number }
  \item{starttime }{character string of date and time at start of fit }
  \item{proctime }{processor time for model fit, in seconds }
 
 The environment variable RCPP_PARALLEL_NUM_THREADS is updated with the 
 value of \code{ncores} if provided.
 
}

\references{

  Gimenez, O., Viallefont, A., Catchpole, E. A., Choquet, R. and Morgan,
  B. J. T. (2004) Methods for investigating parameter redundancy. 
  \emph{Animal Biodiversity and Conservation} \bold{27}, 561--572.

  Huggins, R. M. (1989) On the statistical analysis of capture
  experiments. \emph{Biometrika} \bold{76}, 133--140.

  Pledger, S., Efford, M., Pollock. K., Collazo, J. and Lyons, J. (2009)
  Stopover duration analysis with departure probability dependent on
  unknown time since arrival. In: D. L. Thompson, E. G. Cooch and
  M. J. Conroy (eds) \emph{Modeling Demographic Processes in Marked
  Populations}. Springer. Pp. 349--363.

  Pledger, S., Pollock, K. H. and Norris, J. L. (2010) Open
  capture--recapture models with heterogeneity: II. Jolly-Seber
  model. \emph{Biometrics} \bold{66}, 883--890.

  Pradel, R. (1996) Utilization of capture-mark-recapture for the study
  of recruitment and population growth rate. \emph{Biometrics}
  \bold{52}, 703--709.

  Schwarz, C. J. and Arnason, A. N. (1996) A general methodology for the
  analysis of capture-recapture experiments in open
  populations. \emph{Biometrics} \bold{52}, 860--873.

}

\note{

Different parameterisations lead to different model fits when used with
the default `model' argument in which each real parameter is constrained
to be constant over time.
  
The JSSA implementation uses summation over feasible 'birth' and 'death'
times for each capture history, following Pledger et al. (2010). This
enables finite mixture models for individual capture probability (not
fully tested), flexible handling of additions and losses on capture (aka
removals) (not yet programmed), and ultimately the extension to `unknown
age' as in Pledger et al. (2009).

openCR uses the generalized matrix inverse `ginv' from the MASS
package rather than `solve' from base R, as this seems more robust to
singularities in the Hessian. Also, the default maximization method is `BFGS'
rather than `Newton-Raphson' as BFGS appears more robust in the presence
of redundant parameters.

Earlier versions of \code{\link{openCR.fit}} computed latent class membership 
probabilities for each individual in finite mixture models and saved them in 
component `posterior'. Now see \code{\link{classMembership}} for that functionality.

From 1.5.0 onwards the number of threads uses the environment variable 
RCPP_PARALLEL_NUM_THREADS, as in \code{\link[secr]{secr.fit}}. This may be set once
in a session with \code{secr::setNumThreads}.

The default movement arguments changed in \pkg{openCR} 2.1.1. Now 
\code{kernelradius = 30, sparsekernel = TRUE}.

}

\seealso{
  
  \code{\link{classMembership.openCR}}, 
  \code{\link{derived.openCR}},  
  \code{\link{openCR.design}}, 
  \code{\link{par.openCR.fit}}, 
  \code{\link{predict.openCR}}, 
  \code{\link{summary.openCR}}

}

\examples{

\dontrun{

## CJS default
openCR.fit(ovenCH)

## POPAN Jolly-Seber Schwarz-Arnason, lambda parameterisation
L1 <- openCR.fit(ovenCH, type = 'JSSAl')
predict(L1)

JSSA1 <- openCR.fit(ovenCH, type = 'JSSAf')
JSSA2 <- openCR.fit(ovenCH, type = 'JSSAf', model = list(phi~t))
JSSA3 <- openCR.fit(ovenCH, type = 'JSSAf', model = list(p~t,phi~t))
AIC (JSSA1, JSSA2, JSSA3)
predict(JSSA1)

RMdata <- RMarkInput (join(reduce(ovenCH, by = "all")))
if (require(RMark)) {
    MarkPath <- 'c:/Mark/'
    if (!all (nchar(Sys.which(c('mark.exe', 'mark64.exe', 'mark32.exe'))) < 2)) {
        openCHtest <- process.data(RMdata, model = 'POPAN')
        openCHPOPAN <- mark(data = openCHtest, model = 'POPAN',
            model.parameters = list(p = list(formula = ~1),
            pent = list(formula = ~1),
            Phi = list(formula = ~1)))
        popan.derived(openCHtest, openCHPOPAN)
        cleanup(ask = FALSE)
    } else message ("mark.exe not found")
} else message ("RMark not found")

}

}

\keyword{ model }
