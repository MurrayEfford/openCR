\name{moving.fit}

\alias{moving.fit}
\alias{extractFocal}

\title{Moving Window Functions}

\description{ 

Apply a function to successive multi-session windows from a capthist object. The default function is \code{openCR.fit}, but any function may be used whose first argument accepts a capthist object.

}

\usage{

moving.fit (..., width = 3, centres = NULL, filestem = NULL, 
    trace = FALSE, FUN = openCR.fit)

extractFocal (ocrlist, ...) 

}

\arguments{
  
  \item{...}{named arguments passed to \code{\link{openCR.fit}} (see Details) }
  \item{width}{integer; moving window width (number of primary sessions)}
  \item{centres}{integer; central sessions of windows to consider}
  \item{filestem}{character or NULL; stem used to form filenames for optional  intermediate output}
  \item{trace}{logical; if TRUE a status message is given at each call of FUN}
  \item{FUN}{function to be applied to successive capthist objects}

  \item{ocrlist}{openCRlist object returned by \code{moving.fit} when FUN = openCR.fit}

}

\details{

 \code{moving.fit} applies \code{FUN} to successive multi-session subsets
of the data in the \code{capthist} argument. \code{width} should be an odd integer. 
\code{centres} may be used to restrict the range of windows considered; 
the default is to use all complete windows (width\%/\%2 + 1)...).

If a \code{filestem}
is specified then each result is output to a file that may be loaded with
\code{load}. This is useful if fitting takes a long time and analyses
may be terminated before completion. 

\code{extractFocal} returns the focal-session (central) estimates from a \code{moving.fit} 
with \code{FUN = openCR.fit}. The \dots argument is passed to \code{predict.openCR}; 
it may be used, for example, to choose a different alpha level for confidence intervals.

\code{extractFocal} is untested for complex models (e.g. finite mixtures).

}

\value{

A list in which each component is the output from FUN applied to one subset. The window width is saved as attribute `width'.

}

\seealso{

\code{\link{openCR.fit}}

}

\examples{

## number of individuals detected
moving.fit(capthist = OVpossumCH, FUN = nrow)

\dontrun{

## if package R2ucare installed
if (requireNamespace("R2ucare"))
    moving.fit(capthist = OVpossumCH, FUN = ucare.cjs, width = 5, tests = "overall_CJS")

## using default FUN = openCR.fit

mf1 <- moving.fit(capthist = OVpossumCH, type = 'JSSAfCL', 
     model = list(p~t, phi~t))
lapply(mf1, predict)
extractFocal(mf1)
     
msk <- make.mask(traps(OVpossumCH[[1]]), nx = 32)
mf2 <- moving.fit(capthist = OVpossumCH, mask = msk, type = 'JSSAsecrfCL')
extractFocal(mf2)

}

}