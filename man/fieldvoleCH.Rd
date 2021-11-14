\name{Field vole}
\alias{fieldvoleCH}
\docType{data}

\title{Kielder Field Voles}

\description{

Captures of \emph{Microtus agrestis} on a large grid in a clearcut within Kielder Forest, northern England, June--August 2000 (Ergon and Gardner 2014). Robust-design data from four primary sessions of 3--5 secondary sessions each.

}

\usage{

fieldvoleCH

}

\format{

The format is a multi-session secr capthist object. Attribute `ampm' codes for type of secondary session (am, pm).
  
}

\details{

Ergon and Lambin (2013) provided a robust design dataset from a trapping study on field voles \emph{Microtus agrestis} in a clearcut within Kielder Forest, northern England -- see also Ergon et al. (2011), Ergon and Gardner (2014) and Reich and Gardner (2014). The study aimed to describe sex differences in space-use, survival and dispersal among adult voles. Data were from one trapping grid in summer 2000. 

Trapping was on a rectangular grid of 192 multi-catch (Ugglan Special) traps at 7-metre spacing. Traps were baited with whole barley grains and carrots; voles were marked with individually numbered ear tags. 

Four trapping sessions were conducted at intervals of 21 to 23 days between 10 June and 15 August. Traps were checked at about 12 hour intervals (6 am and 6 pm). 

The attribute `ampm' is a data.frame with a vector of codes, one per secondary session, to separate am and pm trap checks (1 = evening, 2 = morning). The four primary sessions had respectively 3, 5, 4 and 5 trap checks.

Ergon and Gardner (2014) restricted their analysis to adult voles (118 females and 40 males). Histories of five voles (ma193, ma239, ma371, ma143, ma348) were censored part way through the study because they died in traps (T. Ergon pers. comm.).

}

\source{

Data were retrieved from DRYAD (Ergon and Lambin (2013) for \pkg{openCR}. Code for translating the DRYAD ASCII file into a capthist object is given in Examples.

}

\references{

Efford, M. G. (2019) Multi-session models in secr 4.1. \url{https://www.otago.ac.nz/density/pdfs/secr-multisession.pdf}

Ergon, T., Ergon, R., Begon, M., Telfer, S. and Lambin, X. (2011) Delayed density-
dependent onset of spring reproduction in a fluctuating population of field
voles. \emph{Oikos} \bold{120}, 934--940.

Ergon, T. and Gardner, B. (2014) Separating mortality and emigration: modelling space
use, dispersal and survival with robust-design spatial capture--recapture data. 
\emph{Methods in Ecology and Evolution} \bold{5}, 1327--1336.

Ergon, T. and Lambin, X. (2013) Data from: Separating mortality and emigration:
Modelling space use, dispersal and survival with robust-design spatial capture--recapture data.
Dryad Digital Repository. \doi{10.5061/dryad.r17n5}.

Reich, B. J. and Gardner, B. (2014) A spatial capture--recapture model for
territorial species. \emph{Environmetrics} \bold{25}, 630--637.

}

\examples{

summary(fieldvoleCH, terse = TRUE)
m.array(fieldvoleCH)
JS.counts(fieldvoleCH)

attr(fieldvoleCH, 'ampm')

\dontrun{

maleCH <- subset(fieldvoleCH, function(x) covariates(x) == 'M')
fit <- openCR.fit(maleCH)
predict(fit)

# Read data object from DRYAD ASCII file

datadir <- system.file('extdata', package = 'openCR')
EG <- dget(paste0(datadir,'/ergonandgardner2013.rdat'))

# construct capthist object
onesession <- function (sess) {
    mat <- EG$H[,,sess]
    id <- as.numeric(row(mat))
    occ <- as.numeric(col(mat))
    occ[mat<0] <- -occ[mat<0]
    trap <- abs(as.numeric(mat)) 
    matrow <- rownames(mat)
    df <- data.frame(session = rep(sess, length(id)), 
                     ID = matrow[id], 
                     occ = occ, 
                     trapID = trap,
                     sex = c('F','M')[EG$gr],
                     row.names = 1:length(id))
    # retain captures (trap>0)
    df[df$trapID>0, , drop = FALSE]
}
tr <- read.traps(data = data.frame(EG$X), detector = "multi")

# recode matrix as mixture of zeros and trap numbers
EG$H <- EG$H-1

# code censored animals with negative trap number
# two ways to recognise censoring
censoredprimary <- which(EG$K < 4)
censoredsecondary <- which(apply(EG$J,1,function(x) any(x-c(3,5,4,5) < 0)))
censored <- unique(c(censoredprimary, censoredsecondary))
rownames(EG$H)[censored]
# [1] "ma193" "ma239" "ma371" "ma143" "ma348"
censorocc <- apply(EG$H[censored,,], 1, function(x) which.max(cumsum(x)))
censor3 <- ((censorocc-1) \%/\% 5)+1       # session
censor2 <- censorocc - (censor3-1) * 5   # occasion within session
censori <- cbind(censored, censor2, censor3)
EG$H[censori] <- -EG$H[censori] 

lch <- lapply(1:4, onesession)
ch <- make.capthist(do.call(rbind,lch), tr=tr, covnames='sex')

# apply intervals in months
intervals(ch) <-  EG$dt

fieldvoleCH <- ch

# extract time covariate - each secondary session was either am (2) or pm (1)
# EG$tod
# 1 2 3  4  5
# 1 2 1 2 NA NA
# 2 2 1 2  1  1
# 3 2 1 2  1 NA
# 4 2 1 2  1  2
# Note consecutive pm trap checks in session 2
ampm <- split(EG$tod, 1:4)
ampm <- lapply(ampm, na.omit)
attr(fieldvoleCH, 'ampm') <- data.frame(ampm = unlist(ampm))

}

}

\keyword{datasets}
