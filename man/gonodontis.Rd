\name{gonodontisCH}
\alias{gonodontisCH}
\docType{data}

\title{Gonodontis Moths}

\description{
  
Non-spatial open-population capture--recapture data of Bishop et al. (1978) for nonmelanic male \emph{Gonodontis bidentata} at Cressington Park, northwest England.

}

\usage{gonodontisCH}

\format{

  The format is a single-session secr capthist object. As these are
  non-spatial data, the traps attribute is NULL.

}

\details{
  
The data are from a study of the relative fitness of melanic and
nonmelanic morphs of the moth \emph{Gonodontis bidentata} at several
sites in England (Bishop et al. 1978). Crosbie (1979; see also Crosbie
and Manly 1985) selected a subset of the Bishop et al. data (nonmelanic
males from Cressington Park) to demonstrate innovations in Jolly-Seber
modelling, and the same data were used by Link and Barker (2005) and
Schofield and Barker (2008). The present data are those used by Crosbie
(1979) and Link and Barker (2005). 

Male moths were attracted to traps which consisted of a cage containing phermone-producing females surrounded by an enclosure which the males could enter but not leave. New virgin females were usually added every 1 to 4 days. Moths were marked at each capture with a date-specific mark in enamel paint or felt-tip pen on the undersurface of the wing. Thus, although moths at Cressington Park were not marked individually, each moth was a flying bearer of its own capture history.

The data comprise 689 individual capture histories for moths captured at 8 traps operated over 17 days (24 May--10 June 1970). The traps were in a square that appears have been about 40 m on a side. The location of captures is not included in the published data. All captured moths appear to have been marked and released (i.e. there were no removals recorded). All captures on Day 17 were recaptures; it is possible that unmarked moths were not recorded on that day.  

Both Table 1 and Appendix 1 (microfiche) of Bishop et al. (1978) refer to 690 capture histories of nonmelanics at Cressington Park. In the present data there are only 689, and there are other minor discrepancies. Also, Crosbie and Manly (1985: Table 1) refer to 82 unique capture histories (``distinct cmr patterns'') when there are only 81 in the present dataset (note that two moths share 00000000000000011).

}

\source{

Richard Barker provided an electronic copy of the data used by Link and Barker (2005), copied from Crosbie (1979).

}

\references{

Bishop, J. A., Cook, L. M., and Muggleton, J. (1978). The response of two species of moth to
 industrialization in northwest England. II. Relative fitness of morphs
 and population size. \emph{Philosophical Transactions of the Royal Society of London}
 \bold{B281}, 517--540.
 
Crosbie, S. F. (1979) \emph{The mathematical modelling of capture--mark--recapture experiments on animal populations}. Ph.D. Thesis, University of Otago, Dunedin, New Zealand.

Crosbie, S. F. and Manly, B. F. J. (1985) Parsimonious modelling of capture--mark--recapture studies.
\emph{Biometrics} \bold{41}, 385--398.
 
Link, W. A. and Barker, R. J. (2005) Modeling association among demographic parameters
in analysis of open-population capture--recapture data. \emph{Biometrics} \bold{61}, 46--54.

Schofield, M. R. and Barker, R. J. (2008) A unified capture--recapture framework. \emph{Journal of Agricultural
Biological and Environmental Statistics} \bold{13}, 458--477.

}

\examples{

summary(gonodontisCH)
m.array(gonodontisCH)

\dontrun{
# compare default (CJS) estimates from openCR, MARK

fit <- openCR.fit(gonodontisCH)
predict(fit)

if (require(RMark)) {
    MarkPath <- 'c:/Mark/'   # customize as needed
    if (!all (nchar(Sys.which(c('mark.exe','mark64.exe', 'mark32.exe'))) < 2)) {
       mothdf <- RMarkInput(gonodontisCH)
       mark(mothdf)
       cleanup(ask = FALSE)
    } else message ("mark.exe not found")
} else message ("RMark not found")

}

}

\keyword{datasets}
