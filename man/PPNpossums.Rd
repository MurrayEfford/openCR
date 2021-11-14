\name{PPNpossums}
\alias{FebpossumCH}

\docType{data}

\title{

Orongorongo Valley Brushtail Possums

}

\description{

A subset of brushtail possum (\emph{Trichosurus vulpecula}) data from the Orongorongo Valley live-trapping study of Efford (1998) and Efford and Cowan (2005) that was used by Pledger, Pollock and Norris (2003, 2010). The \code{OVpossumCH} dataset in \pkg{secr} is a different selection of data from the same study. Consult ?OVpossumCH for more detail.

The data comprise captures in February of each year from 1980 to 1988. 
}

\usage{

FebpossumCH

}

\format{

The format is a 9-session \pkg{secr} capthist object. Capture locations are not included.

}

\details{

The data are captures of 448 animals (175 females and 273 males) over 9 trapping sessions comprising 4--10 occasions each. All were independent of their mothers, but age was not otherwise distinguished. The individual covariate \code{sex} takes values `F' or `M'.

Pledger, Pollock and Norris (2010) fitted 2-class finite mixture models for capture probability p and apparent survival phi, with or without allowance for temporal (between year) variation, using captures from only the first day of each trapping session. The first-day data relate to 270 individuals (115 females and 155 males). 

}

\source{

M. Efford unpubl. See Efford and Cowan (2004) for acknowledgements.

}

\references{

Efford, M. G. (1998) Demographic consequences of sex-biased dispersal in
     a population of brushtail possums. \emph{Journal of Animal Ecology}
     \bold{67}, 503--517.

Efford, M. G. and Cowan, P. E. (2004) Long-term population trend of
     \emph{Trichosurus vulpecula} in the Orongorongo Valley, New
     Zealand. In: \emph{The Biology of Australian Possums and
     Gliders}. Edited by R. L. Goldingay and S. M. Jackson. Surrey
     Beatty & Sons, Chipping Norton. Pp. 471--483.
     
Pledger, S., Pollock, K. H. and Norris, J. L. (2010) Open
  capture--recapture models with heterogeneity: II. Jolly--Seber
  model. \emph{Biometrics} \bold{66}, 883--890.
  
}

\examples{

summary(FebpossumCH) 
m.array(FebpossumCH)
JS.counts(FebpossumCH)

FebD1CH <- subset(FebpossumCH, occasion = 1)

\dontrun{

# reading the text file 'poss8088.data'

datadir <- system.file('extdata', package = 'openCR')
poss8088df <- read.table (paste0(datadir,'/poss8088.data'), header = TRUE)
capt <- poss8088df[,c('session','id','day','day','sex')]

# duplication of day is a trick to get a dummy trapID column in the right place
# this is needed because make.capthist does not have nonspatial option
capt$day.1[] <- 1  

# keep only February samples
capt <- capt[capt$session \%\% 3 == 1,]

# build nonspatial secr capthist object using dummy trapping grid
FebpossumCH <- make.capthist(capt, make.grid(1,2,ID='numx'))
# discard dummy traps objects
for (i in 1:9) attr(FebpossumCH[[i]], 'traps') <- NULL
names(FebpossumCH) <- 1980:1988 
sessionlabels(FebpossumCH) <- 1980:1988

}

}
\keyword{datasets}
