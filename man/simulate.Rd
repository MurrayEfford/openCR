\name{simulation}
\alias{sim.nonspatial}
\alias{runsim.nonspatial}
\alias{runsim.spatial}
\alias{runsim.RMark}
\alias{sumsims}

\title{Simulate Capture Histories}

\description{

Generate non-spatial or spatial open-population data and fit models.

}

\usage{

sim.nonspatial (N, turnover = list(), p, nsessions, noccasions = 1, intervals = NULL, 
    recapfactor = 1, seed = NULL, savepopn = FALSE, ...)
    
runsim.nonspatial (nrepl = 100, seed = NULL, ncores = NULL, fitargs = list(), 
    extractfn = predict, ...)

runsim.spatial (nrepl = 100, seed = NULL, ncores = NULL, popargs = list(), 
    detargs = list(), fitargs = list(), extractfn = predict, intervals = NULL)

sumsims (sims, parm = 'phi', session = 1, dropifnoSE = TRUE, svtol = NULL, 
    maxcode = 3, true = NULL)

runsim.RMark (nrepl = 100, model = "CJS", model.parameters = NULL, extractfn,
    seed = NULL, ...)

}
\arguments{

  \item{N}{integer population size}
  \item{turnover}{list as described for \link[secr]{turnover}}
  \item{p}{numeric detection probability}
  \item{nsessions}{number of primary sessions}
  \item{noccasions}{number of secondary sessions per primary session}
  \item{intervals}{intervals between primary sessions (see Details)}
\item{recapfactor}{numeric multiplier for capture probability after
  first capture}
  \item{seed}{random number seed see \link[secr]{random numbers}}
  \item{savepopn}{logical; if TRUE the generated population is saved as an attribute of the capthist object}
  \item{\dots}{other arguments passed to \code{\link[secr]{sim.popn}} (sim.nonspatial) or \code{\link{sim.nonspatial}} (runsims)}
  
  \item{nrepl}{number of replicates}
  \item{ncores}{ integer number of cores to be used for parallel processing (see Details)}
  \item{popargs}{list of arguments for sim.popn}
  \item{detargs}{list of arguments for sim.capthist}
  \item{fitargs}{list of arguments for openCR.fit}
  \item{extractfn}{function applied to each fitted openCR model}
  
  \item{sims}{list output from \code{runsim.nonspatial} or \code{runsim.spatial}}
  \item{parm}{character name of parameter to summarise}
  \item{session}{integer vector of session numbers to summarise}
  \item{dropifnoSE}{logical; if TRUE then replicates are omitted when SE missing for parm}
  \item{svtol}{numeric; minimum singular value (eigenvalue) considered non-zero}
  \item{maxcode}{integer; maximum accepted value of convergence code}
  \item{true}{true value of requested parm in given session}
  
  \item{model}{character; RMark model type }
  \item{model.parameters}{list with RMark model specification (see \code{?mark})}

}

\details{

For \code{sim.nonspatial} -- If \code{intervals} is specified then the number of primary and secondary sessions is inferred from \code{intervals} and \code{nsessions} and \code{noccasions} are ignored. If \code{N} and \code{p} are vectors of length 2 then subpopulations of the given initial size are sampled with the differing capture probabilities and the resulting capture histories are combined.

\code{runsim.spatial} is a relatively simple wrapper for \code{\link[secr]{sim.popn}}, \code{\link[secr]{sim.capthist}}, and \code{\link{openCR.fit}}. Some arguments are set automatically: the \code{sim.capthist} argument 'renumber' is always FALSE; argument 'seed' is ignored within 'popargs' and 'detargs'; if no 'traps' argument is provided in 'detargs' then 'core' from 'popargs' will be used; detargs$popn and fitargs$capthist are derived from the preceding step. The 'type' specified in fitargs may refer to a non-spatial or spatial open-population model ('CJS', 'JSSAsecrfCL' etc.). If the \code{intervals} argument is specified it is used to set the intervals attribute of the simulated capthist object; turnover parameters in \code{sim.popn} are not scaled by \code{intervals}.

Control of parallel processing changed in \pkg{openCR} 1.5.0 to conform to \pkg{secr}. In \code{runsim.nonspatial} and \code{runsim.spatial}, if \code{ncores} is NULL (the default) then the number of cores used for multithreading by \code{openCR.fit} is controlled by the environment variable RCPP_PARALLEL_NUM_THREADS. Use the secr function \code{setNumThreads} to set RCPP_PARALLEL_NUM_THREADS to a value greater than the default (2, from \pkg{openCR} 1.5 onwards).

Otherwise, (\code{ncores} specified in runsim.nonspatial or runsim.spatial) 'ncores' is set to 1 for each replicate and the replicates are split across the specified number of cores.

\code{sumsims} assumes output from \code{runsim.nonspatial} and \code{runsim.spatial} with `extractfn = predict' or `extractfn = summary'. Missing SE usually reflects non-identifiability of a parameter or failure of maximisation, so these replicates are dropped by default. If \code{svtol} is specified then the rank of the Hessian is determined by counting eigenvalues that exceed svtol, and replicates are dropped if the rank is less than the number of beta parameters. A value of 1e-5 is suggested for svtol in \code{\link{AIC.openCR}}, but smaller values may be appropriate for larger models (MARK has its own algorithm for this threshold).

Replicates are also dropped by \code{sumsims} if the convergence code exceeds 'maxcode'. The maximisation functions \code{\link{nlm}} (used for method = 'Newton-Raphson', the default), and \code{\link{optim}} (all other methods) return different convergence codes; their help pages should be consulted. The default is to accept code = 3 from \code{nlm}, although the status of such maximisations is ambiguous.

}

\value{

\code{sim.nonspatial} --

A capthist object representing an open-population sample

\code{runsim.nonspatial} and \code{runsim.spatial} --

List with one component (output from extractfn) for each replicate. Each component also has attributes 'eigH' and 'fit' taken from the output of \code{openCR.fit}. See Examples to extract convergence codes from 'fit' attribute.

\code{sumsims} --

Data.frame with rows `estimate', `SE.estimate', `lcl', `ucl', `RSE', `CI.length' and columns for median, mean, SD and n. If `true' is specified there are additional rows are `Bias' and `RB', and columns for `rRMSE' and `COV'.

}

\seealso{

\code{\link[secr]{sim.popn}}, \code{\link[secr]{sim.capthist}}

}

\examples{

\dontrun{

cores <- 2   # for CRAN check; increase as available

ch <- sim.nonspatial(100, list(phi = 0.7, lambda = 1.1), p = 0.3, nsessions = 8, noccasions=2)
openCR.fit(ch, type = 'CJS')

turnover <- list(phi = 0.85, lambda = 1.0, recrmodel = 'constantN')
set.seed(123)

## using type = 'JSSAlCL' and extractfn = predict

fitarg <- list(type = 'JSSAlCL', model = list(p~t, phi~t, lambda~t))
out <- runsim.nonspatial(nrepl = 100, N = 100, ncores = cores, turnover = turnover, 
   p = 0.2, recapfactor = 4, nsessions = 10, noccasions = 1, fitargs = fitarg)
sumsims(out, 'lambda', 1:10)

## using type = 'Pradelg' and extractfn = derived
## homogeneous p
fitarg <- list(type = 'Pradelg', model = list(p~t, phi~t, gamma~t))
outg <- runsim.nonspatial(nrepl = 100, N = 100, ncores = cores, turnover = turnover, 
    p = 0.2, recapfactor = 4, nsessions = 10, noccasions = 1, 
    fitargs = fitarg, extractfn = derived)
apply(sapply(outg, function(x) x$estimates$lambda),1,mean)

turnover <- list(phi = 0.85, lambda = 1.0, recrmodel = 'discrete')

## 2-class mixture for p
outg2 <- runsim.nonspatial(nrepl = 100, N = c(50,50), ncores = cores, turnover = turnover, 
    p = c(0.3,0.9), recapfactor = 1, nsessions = 10, noccasions = 1, 
    fitargs = fitarg, extractfn = derived)
outg3 <- runsim.nonspatial(nrepl = 100, N = c(50,50), ncores = cores, turnover = turnover, 
    p = c(0.3,0.3), recapfactor = 1, nsessions = 10, noccasions = 1, 
    fitargs = fitarg, extractfn = derived)
apply(sapply(outg2, function(x) x$estimates$lambda),1,mean)

plot(2:10, apply(sapply(outg2, function(x) x$estimates$lambda),1,mean)[-1],
    type='o', xlim = c(1,10), ylim = c(0.9,1.1))

## RMark

extfn <- function(x) x$results$real$estimate[3:11]
MarkPath <- 'c:/mark'  # customise as needed
turnover <- list(phi = 0.85, lambda = 1.0, recrmodel = 'discrete')
outrm <- runsim.RMark (nrepl = 100, model = 'Pradlambda', extractfn = extfn, 
                       model.parameters = list(Lambda=list(formula=~time)),
                       N = c(200,200), turnover = turnover, p = c(0.3,0.9),
                       recapfactor = 1, nsessions = 10, noccasions = 1)
extout <- apply(do.call(rbind, outrm),1,mean)

## Spatial

grid <- make.grid()
msk <- make.mask(grid, type = 'trapbuffer', nx = 32)
turnover <- list(phi = 0.8, lambda = 1)
poparg <- list(D = 10, core = grid, buffer = 100, Ndist = 'fixed', nsessions = 6, 
    details = turnover)
detarg <- list(noccasions = 5, detectfn = 'HHN', detectpar = list(lambda0 = 0.5, sigma = 20))
fitarg <- list(type = 'JSSAsecrfCL', mask = msk, model = list(phi~1, f~1))
sims <- runsim.spatial (nrepl = 7, ncores = cores, pop = poparg, det = detarg, fit = fitarg)
sumsims(sims)

## extract the convergence code from nlm for each replicate in preceding simulation
sapply(lapply(sims, attr, 'fit'), '[[', 'code')
## if method != 'Newton-Raphson then optim is used and the code is named 'convergence'
# sapply(lapply(sims, attr, 'fit'), '[[', 'convergence')

}

}

\keyword{ datagen }