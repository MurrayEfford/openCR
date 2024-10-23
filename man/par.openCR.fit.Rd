\name{par.openCR.fit}
\alias{par.openCR.fit}

\title{Fit Multiple openCR Models}

\description{
  
  This function is a wrapper for \code{\link{openCR.fit}}. 
  
}

\usage{

par.openCR.fit (arglist, ncores = 1, seed = 123, trace = FALSE, logfile = NULL, 
    prefix = "")

}

\arguments{
  
  \item{arglist}{list of argument lists for \code{secr.fit} or a character
    vector naming such lists}
  
  \item{ncores}{ integer number of cores used by parallel::makeClusters() }

  \item{seed}{integer pseudorandom number seed}

  \item{trace}{logical; if TRUE intermediate output may be logged}
  
  \item{logfile}{character name of file to log progress reports}
  
  \item{prefix}{character prefix for names of output}
  
}

\details{

In openCR >= 1.5.0, setting ncores > 1 is deprecated and triggers a warning: multithreading makes it faster to set ncores = 1 in par.openCR.fit.

\code{trace} overrides any settings in \code{arglist}.
  
It is convenient to provide the names of the capthist and mask arguments
  in each component of arglist as character values (i.e. in quotes); objects thus
  named are exported from the workspace to each worker process (see Examples).

Using \code{ncores}>1 is obsolete under the multithreading regime in \pkg{openCR} >= 1.5.0. 
It is usually slower than \code{ncores} = 1. If used it has these effects:

-- worker processes are generated using the \pkg{parallel} package,

-- one model is fitted on each worker, and

-- if no logfile name is provided then a temporary file name will be generated in tempdir().

}

\value{

For \code{par.openCR.fit} - openCRlist of model fits (see
\code{\link{openCR.fit}} and \code{\link{openCRlist}}). Names are created by prefixing \code{prefix} to the
names of \code{argslist}. If \code{trace} is TRUE then the total
execution time and finish time are displayed.

}

\seealso{
  \code{\link{openCR.fit}},
  \link[secr]{Parallel},
  \code{\link{make.table}},
  \code{\link{openCRlist}}
}

\note{

Any attempt in \code{arglist} to set \code{ncores > 1} for a particular 
openCR fit was ignored in \pkg{openCR} < 1.5.0. Now it is allowed.

}

\examples{

\dontrun{

m1 <- list(capthist = ovenCH, model = list(p~1, phi~1)) 
m2 <- list(capthist = ovenCH, model = list(p~session, phi~1))
m3 <- list(capthist = ovenCH, model = list(p~session, phi~session) )
setNumThreads(7)  # on quadcore Windows PC
fits <- par.openCR.fit (c('m1','m2','m3'), ncores = 1)
AIC(fits)

}
}
\keyword{ model }

