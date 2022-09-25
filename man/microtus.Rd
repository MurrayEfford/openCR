\name{Microtus}
\alias{microtusCH}
\alias{microtusFCH}
\alias{microtusMCH}
\alias{microtusFMCH}
\alias{microtusRDCH}
\docType{data}

\title{Patuxent Meadow Voles}

\description{

Captures of \emph{Microtus pennsylvanicus} at Patuxent Wildlife Research Center, Laurel, Maryland, June--December 1981. Collapsed (primary session only) data for adult males and adult females, and full robust-design data for adult males. Nichols et al. (1984) described the field methods and analysed a superset of the present data.

}

\usage{

microtusCH
microtusFCH
microtusMCH
microtusFMCH
microtusRDCH

}

\format{

  The format is a single-session secr capthist object. As these are
  non-spatial data, the traps attribute is NULL.
  
}

\details{

Voles were caught in live traps on a 10 x 10 grid with traps 7.6 m apart. Traps were baited with corn. Traps were set in the evening, checked the following morning, and locked open during the day. Voles were ear-tagged with individually numbered fingerling tags. The locations of captures were not included in the published data. 

Data collection followed Pollock's robust design with five consecutive days of trapping each month for six months (27 June 1981--8 December 1981). The data are for "adult" animals only, defined as those weighing at least 22g. Low capture numbers on the last two days of the second primary session (occasions 9 and 10) are due to a raccoon interfering with traps (Nichols et al. 1984). Six adult female voles and ten adult male voles were not released; their final captures are coded as -1 in the respective capthist objects.

\code{microtusRDCH} is the full robust-design dataset for adult males ((Williams et al. 2002 Table 19.1). 

\code{microtusFCH} and \code{microtusMCH} are the collapsed datasets (binary at the level of primary session) for adult females and adult males from Williams et al. (2002 Table 17.5); \code{microtusFMCH} combines them and includes the covariate `sex'.

\code{microtusCH} is a combined-sex version of the data with different lineage (see below).

The `intervals' attribute was assigned for \code{microtusRDCH} to distinguish primary sesssions (interval 1 between prmary sessions; interval 0 for consecutive secondary sessions within a primary session). True intervals (start of one primary session to start of next) were 35, 28, 35, 28 and 34 days. See Examples to add these manually.

Williams, Nichols and Conroy (2002) presented several analyses of these data. 

Program JOLLY (Hines 1988, Pollock et al. 1990) included a combined-sex version of the primary-session data that was used by Pollock et al. (1985) and Pollock et al. (1990)*. The numbers of voles released each month in the JOLLY dataset JLYEXMPL differ by 0--3 from the sum of the male and female data from Williams et al. (2002) (see Examples). Some discrepancies may have been due to voles for which sex was not recorded. The JOLLY version matches Table 1 of Nichols et al. (1984). The JOLLY version is distributed here as the object \code{microtusCH}.

Differing selections of data from the Patuxent study were analysed by Nichols et al. (1992) and Bonner and Schwarz (2006).

* There is a typographic error in Table 4.7 of Pollock et al. (1990): \eqn{r_i} for the first period should be 89.

}

\source{

\tabular{llll}{
Object \tab  \tab Source \cr
\code{microtusCH} \tab  \tab Text file JLYEXMPL distributed with Program JOLLY (Hines 1988; see also Examples) \cr
\code{microtusFCH} \tab  \tab Table 17.5 in Williams, Nichols and Conroy (2002) \cr
\code{microtusMCH} \tab  \tab Table 17.5 in Williams, Nichols and Conroy (2002) \cr
\code{microtusFMCH} \tab  \tab Table 17.5 in Williams, Nichols and Conroy (2002) \cr
\code{microtusRDCH} \tab  \tab Table 19.1 in Williams, Nichols and Conroy (2002) provided as text file by Jim Hines\cr
}

}

\references{

Bonner, S. J. and Schwarz, C. J. (2006) An extension of the Cormack--Jolly--Seber model for continuous
covariates with application to \emph{Microtus pennsylvanicus}. \emph{Biometrics} \bold{62}, 142--149.

Hines, J. E. (1988) Program "JOLLY". Patuxent Wildlife Research Center. \url{https://www.mbr-pwrc.usgs.gov/software/jolly.html}

Nichols, J. D., Pollock, K. H., Hines, J. E. (1984) The use of a robust capture-recapture design in small mammal population studies: a field example with \emph{Microtus pennsylvanicus}. \emph{Acta Theriologica} \bold{29}, 357--365.

Nichols, J. D., Sauer, J. R., Pollock, K. H., and Hestbeck, J. B. (1992) Estimating transition probabilities for stage-based population projection matrices using capture--recapture
data. \emph{Ecology} \bold{73}, 306--312.

Pollock, K. H., Hines, J. E. and Nichols, J. D. (1985) Goodness-of-fit tests for open capture--recapture models. 
\emph{Biometrics} \bold{41}, 399--410.

Pollock, K. H., Nichols, J. D., Brownie, C. and Hines, J. E. (1990) Statistical inference for capture--recapture experiments. \emph{Wildlife Monographs} \bold{107}. 97pp. 

Williams, B. K., Nichols, J. D. and Conroy, M. J. (2002) \emph{Analysis and management of animal populations.} Academic Press.

}

\examples{

# cf Williams, Nichols and Conroy Table 17.6
m.array(microtusFCH)
m.array(microtusMCH)

\dontrun{

# cf Williams, Nichols and Conroy Fig. 17.2
fitfm <- openCR.fit(microtusFMCH, model = list(p~1, phi ~ session + sex))
maledat <- expand.grid(sex = factor('M', levels = c('F','M')), session = factor(1:6))
plot(fitfm, ylim=c(0,1), type = 'o')
plot(fitfm, newdata = maledat, add = TRUE, xoffset = 0.1, pch = 16, type = 'o')

# adjusting for variable interval
intervals(microtusCH) <-  c(35,28,35,28,34) / 30 
intervals(microtusRDCH)[intervals(microtusRDCH)>0] <- c(35,28,35,28,34) / 30

# The text file JLYEXMPL distributed with JOLLY is in the extdata folder of the R package
# The microtusCH object may be rebuilt as follows
datadir <- system.file('extdata', package = 'openCR')
JLYdf <- read.table(paste0(datadir,'/JLYEXMPL'), skip = 3, 
                    colClasses = c('character','numeric'))
names(JLYdf) <- c('ch', 'freq')
JLYdf$freq[grepl('2', JLYdf$ch)] <- -JLYdf$freq[grepl('2', JLYdf$ch)]
JLYdf$ch <- gsub ('2','1', JLYdf$ch)
microtusCH <- unRMarkInput(JLYdf)

# Compare to combined-sex data from Williams et al. Table 17.5
JS.counts(microtusCH) - JS.counts(microtusFMCH)
}

}

\keyword{datasets}
