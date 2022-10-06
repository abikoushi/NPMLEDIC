\name{eccdf_DICRT}
\alias{eccdf_DICRT}
\title{
Estimate emperical survival function from doubly interval censored with right truncated data
}
\description{
Estimate emperical survival function from doubly interval censored with right truncated data
}
\usage{
eccdf_DICRT(EL, ER, SL, SR, tmax=Inf, breaks=NULL, midp=0.5, iter=1000L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EL}{
explosure left
}
  \item{ER}{
explosure right
}
  \item{SL}{
symptom onset left
}
  \item{SR}{
symptom onset right
}
}
\details{
Estimate emperical survival function from doubly interval censored with right truncated data
}
\value{
  \item{time}{time}
  \item{prob}{probability}
}
\examples{
x <- rweibull(100,2,7)
ew <- rexp(100)
sw <- rexp(100)
dat <- simDIC(x,ew,sw)
out <- eccdf_DICRT(EL=dat$E_L, ER=dat$E_R, SL=dat$S_L, SR=dat$S_R)
plot(out$time, out$prob, type="s", ylab="eccdf", xlab="x")
curve(pweibull(x,2,7,lower.tail=FALSE), add=TRUE, col="royalblue")
}