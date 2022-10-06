\name{simDIC}
\alias{simDIC}
\title{
Simulate doubly interval censored
}
\description{
Simulate doubly interval censored
}
\usage{
simDIC(inc, Ew, Sw)
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