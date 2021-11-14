\name{read.inp}
\alias{read.inp}

\title{ Import Data from RMark Input Format }

\description{
  
\code{read.inp} forms a capthist object from a MARK input (.inp) file.

}

\usage{

read.inp(filename, ngroups = 1, grouplabel = 'group', grouplevels = NULL, 
    covnames = NULL, skip = 0)

}

\arguments{

  \item{filename}{character file name including `.inp'.}
  \item{ngroups}{integer number of group columns in input }
  \item{grouplabel}{character}
  \item{grouplevels}{vector with length equal to number of groups}
  \item{covnames}{character vector of additional covariates names, one per covariate column }
  \item{skip}{integer number of lines to skip at start of file }

}
\details{

Comments bracketed with `/*` and `*/' will be removed automatically.

If \code{grouplevels} is specified then \code{ngroups} is taken from the number of levels (\code{ngroups} is overridden). An individual covariate is output, named according to \code{grouplabel}. The order of levels in \code{grouplevels} should match the order of the group frequency columns in the input. This also determines the ordering of levels in the resulting covariate.

}

\value{

  A single-session capthist object with no traps attribute.

}

\seealso{

\code{\link{RMarkInput}}, \code{\link{unRMarkInput}}

}

\examples{

datadir <- system.file('extdata', package = 'openCR')
dipperCH <- read.inp(paste0(datadir, '/ed.inp'), ngroups = 2)
summary(dipperCH)

}

\keyword{ manip }
