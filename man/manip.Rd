\name{miscellaneous}
\alias{primarysessions}
\alias{secondarysessions}
%\alias{timevaryingcov}
%\alias{timevaryingcov<-}

\title{Data Manipulation}

\description{

Miscellaneous functions

}

\usage{

primarysessions(intervals)
secondarysessions(intervals)
%timevaryingcov(object, ...)
% timevaryingcov(object) <- value

}

\arguments{

  \item{intervals}{numeric vector of intervals for time between
    secondary sessions a of robust design}
%  \item{object}{single-session capthist object}
%  \item{value}{a list of named vectors}
%  \item{...}{other arguments (not used)}
  
}

\details{

These functions are used internally.

}

\value{

primarysessions --

Integer vector with the number of the primary session to which each secondary session belongs.

secondarysessions --

Integer vector with secondary sessions numbered sequentially within primary sessions.

}

\examples{

int <- intervals(join(ovenCH))
primary <- primarysessions(int)
primary

# number of secondary sessions per primary
table(primary) 

# secondary session numbers
secondarysessions(int)

}

\keyword{ manip }