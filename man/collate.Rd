\name{collate}
\alias{collate}
\alias{collate.openCR}
\alias{collate.openCRlist}

\title{Array of Parameter Estimates}

\description{
  
Estimates from one or more openCR models are formed into an array.

}

\usage{

\method{collate}{openCR}(object, ..., realnames = NULL, betanames = NULL, 
    newdata = NULL, alpha = 0.05, perm = 1:4, fields = 1:4)

\method{collate}{openCRlist}(object, ..., realnames = NULL, betanames = NULL, 
    newdata = NULL, alpha = 0.05, perm = 1:4, fields = 1:4)
    
}

\arguments{

  \item{object}{\code{openCR} or \code{openCRlist} objects}
  \item{\dots}{ other \code{openCR} objects}
  \item{realnames}{ character vector of real parameter names }
  \item{betanames}{ character vector of beta parameter names }
  \item{newdata}{ optional dataframe of values at which to evaluate models }
  \item{alpha}{ alpha level for confidence intervals }
  \item{perm}{ permutation of dimensions in output }
  \item{fields}{vector to restrict summary fields in output }

}

\details{

\code{collate} extracts parameter estimates from a set of fitted openCR
model objects. \code{fields} may be used to select a subset of summary
fields ("estimate","SE.estimate","lcl","ucl") by name or number. 

}

\value{

A 4-dimensional array of model-specific parameter
estimates. By default, the dimensions correspond respectively to rows in
\code{newdata} (usually sessions), models, statistic fields (estimate, SE.estimate, lcl,
ucl), and parameters ("phi", "sigma" etc.). For particular comparisons it often helps
to reorder the dimensions with the \code{perm} argument.

}

\seealso{ 

\code{\link{modelAverage.openCR}}, 
\code{\link{make.table}}

}

