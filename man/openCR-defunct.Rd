\name{openCR-defunct}
\alias{openCR-defunct}

\title{Defunct Functions in Package \pkg{openCR}}

\description{

These functions are no longer available in \pkg{openCR}.

}

\usage{

# Defunct in 2.2.3

ucare.cjs()

# Defunct in 2.1.0

openCR.make.newdata()

# Defunct in 2.0.0

plotKernel()

}

\details{

\code{ucare.cjs()} was removed from \pkg{openCR} 2.2.4 because \pkg{R2ucare} 
upon which it depends is not currently available from CRAN.

Internal function \code{openCR.make.newdata} was replaced with a method for 
the openCR class of the generic \code{makeNewData}.

\code{plotKernel} was replaced with a plot method for the kernel class.

}

\seealso{

\code{\link{openCR-deprecated}}

}
