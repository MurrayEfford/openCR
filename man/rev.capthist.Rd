\name{rev.capthist}
\alias{rev.capthist}

\title{Reverse Primary Sessions}

\description{

The \code{rev} method for capthist objects reverses the order of the primary sessions while retaining the order of secondary sessions within each primary session.

}

\usage{

\method{rev}{capthist}(x)

}

\arguments{

  \item{x}{multi-session capthist object from secr}
  
}

\details{

rev() is used to demonstrate 'reversed time' analyses (Nichols 2016)
in which seniority (gamma) is estimated as reversed-time survival (phi)
The approach is numerically equivalent to direct modelling of seniority (see Examples). 
Direct modelling allows more control and is more intuitive.

If \code{x} is not overtly multi-session and has no intervals attribute then each occasion is treated as a primary session.

}

\value{

Capthist object with same observations as input, but re-ordered. 
The order of attributes \code{sessionlabels} and \code{intervals} is also reversed.
A default intervals attribute is added if the input lacks one.

}

\references{

Nichols, J. D. (2016) And the first one now will later be last:
time-reversal in Cormack--Jolly--Seber Models. \emph{Statistical Science}
\bold{31,} 175--190.

}

\examples{

summary(rev(ovenCH), terse = TRUE)

# These three models give the same result for gamma except for
# gamma(1982) which is confounded with p and not separately estimable:

\dontrun{

dipperPradel <- openCR.fit(dipperCH, type = "Pradelg", model = list(p~t, phi~t, gamma~t))
revdipper <- openCR.fit(rev(dipperCH), model=list(p~t, phi~t))
dipperJSSA <- openCR.fit(dipperCH, type='JSSAgCL', model=list(p~t, phi~t, gamma~t))

predict(dipperPradel)$gamma
predict(revdipper)$phi
predict(dipperJSSA)$gamma

}

}

\keyword{ manip }