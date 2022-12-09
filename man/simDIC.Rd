\name{simDIC}
\alias{simDIC}
\title{
Simulate doubly interval censored data
}
\description{
Simulate doubly interval censored data
}
\usage{
simDIC(inc, Ew, Sw)
}

\arguments{
  \item{inc}{
  incubation period (the target of estimation)
}
\item{at}{
 time-points of infection (if \code{NULL} draw from \code{rexp})
}
\item{WE}{
 width of the exposure intervals (if \code{NULL} draw from \code{rexp})
}
\item{WS}{
 width of the sympton onset intervals (if \code{NULL} draw from \code{rexp})
}

}
\details{
Simulate doubly interval censored data
}

\examples{
x <- rweibull(10,2,7)
simDIC(x)
}