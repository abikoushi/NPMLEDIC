\name{simDIC}
\alias{simDIC}
\alias{simIC}
\title{
Simulate doubly interval censored data
}
\description{
Simulate doubly interval censored data
}
\usage{
simDIC(inc, Ew, Sw)

simIC(inc, Ew)
}

\arguments{
  \item{inc}{
  incubation period (target of estimation)
}
\item{Ew}{
 width of the exposure intervals
}
\item{Sw}{
 width of the sympton onset intervals
}
}
\details{
simulate
}

\examples{
x <- rweibull(10,2,7)
ew <- rexp(10)
sw <- rexp(10)
simDIC(x,ew,sw)
}