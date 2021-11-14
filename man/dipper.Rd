\name{dipperCH}
\alias{dipperCH}

\docType{data}

\title{Dippers}

\description{

Lebreton et al. (1992) demonstrated Cormack-Jolly-Seber methods with a dataset on European Dipper (*Cinclus cinclus*) collected by Marzolin (1988) and the data have been much used since then. Dippers were captured annually over 1981--1987. We use the version included in the RMark package (Laake 2013).

}

\usage{

dipperCH

}

\format{

 The format is a single-session secr capthist object. As these are
  non-spatial data, the traps attribute is NULL.
 
}

\details{

Dippers were sampled in 1981--1987.

}

\source{

MARK example dataset `ed.inp'. Also RMark (Laake 2013). See Examples.

}

\seealso{

\code{\link{read.inp}}

}

\references{

Laake, J. L. (2013). \emph{RMark: An R Interface for Analysis of Capture--Recapture Data with MARK}.
AFSC Processed Report 2013-01, 25p. Alaska Fisheries Science Center, NOAA, 
National Marine Fisheries Service, 7600 Sand Point Way NE, Seattle WA 98115.

Lebreton, J.-D., Burnham, K. P., Clobert, J., and Anderson,
D. R. (1992) Modeling survival and testing biological hypotheses using
marked animals: a unified approach with case studies. \emph{Ecological
Monographs} \bold{62}, 67--118.

Marzolin, G. (1988) Polygynie du Cincle plongeur (*Cinclus cinclus*) dans 
les c?tes de Lorraine. \emph{L'Oiseau et la Revue Francaise
  d'Ornithologie} \bold{58}, 277--286.

}

\examples{

m.array(dipperCH)

\dontrun{

# From file 'ed.inp' in MARK input format
datadir <- system.file('extdata', package = 'openCR')
dipperCH <- read.inp(paste0(datadir, '/ed.inp'), grouplabel='sex',
    grouplevels = c('Male','Female'))
intervals(dipperCH) <- rep(1,6)    
sessionlabels(dipperCH) <- 1981:1987   # labels only

# or extracted from the RMark package with this code
if (require(RMark)) {
    if (all (nchar(Sys.which(c('mark.exe','mark64.exe', 'mark32.exe'))) < 2))
        stop ("MARK executable not found; set e.g. MarkPath <- 'c:/Mark/'")
    data(dipper)                       # retrieve dataframe of dipper capture histories
    dipperCH2 <- unRMarkInput(dipper)  # convert to secr capthist object
    intervals(dipperCH2) <- rep(1,6)    
    sessionlabels(dipperCH2) <- 1981:1987   # labels only
} else message ("RMark not found")

# The objects dipperCH and dipperCH2 differ in the order of factor levels for 'sex'
}

}

\keyword{datasets}
