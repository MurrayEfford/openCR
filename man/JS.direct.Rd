\name{JS.direct}
\alias{JS.direct}

\title{Jolly--Seber Estimates}

\description{

Non-spatial open-population estimates using the conventional closed-form Jolly--Seber estimators (Pollock et al. 1990).

}

\usage{

JS.direct(object)

}

\arguments{

  \item{object}{secr capthist object or similar}

}

\details{

Estimates are the session-specific Jolly-Seber estimates with no constraints.

The reported SE of births (B) differ slightly from those in Pollock et al. (1990), and may be in error.

}

\value{

A dataframe in which the first 5 columns are summary statistics (counts from \code{\link{JS.counts}}) and the remaining columns are estimates:  

\item{p}{capture probability}
\item{N}{population size}
\item{phi}{probability of survival to next sample time}
\item{B}{number of recruits at next sample time}

Standard errors are in fields prefixed `se'; for N and B these include only sampling variation and omit population stochasticity. The covariance of successive phi-hat is in the field `covphi'.

}

\references{

Pollock, K. H., Nichols, J. D., Brownie, C. and Hines, J. E. (1990) Statistical inference for capture--recapture experiments. \emph{Wildlife Monographs} \bold{107}. 97pp. 

}

\seealso{

\code{\link{JS.counts}}

}

\examples{

# cf Pollock et al. (1990) Table 4.8
JS.direct(microtusCH)

}

