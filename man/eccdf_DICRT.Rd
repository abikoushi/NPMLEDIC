\name{eccdf_DICRT}
\alias{eccdf_DICRT}
\alias{eccdf_ICRT}
\title{
Estimate emperical survival function from doubly interval censored with right truncated data
}
\description{
Estimate emperical survival function from doubly interval censored with right truncated data
}
\usage{
eccdf_DICRT(EL, ER, SL, SR, tmax=Inf, breaks=NULL, midp=0.5, iter=1000L)

eccdf_ICRT(EL, ER, S, tmax=Inf, breaks=NULL, iter=1000L)
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
  \item{S}{
symptom onset
}
\item{tmax}{
upper truncation point
}
\item{breaks}{
breaks of histrogram
}
\item{iter}{
nomber of iteration
}
}
\details{
Estimate emperical survival function from doubly interval censored with right truncated data. 
}
\examples{
set.seed(123456)
n <- 100L
x <- rweibull(n,1.5,5)
ew <- rexp(n)
sw <- rexp(n)
dat <- simDIC(x, ew, sw)
out <- eccdf_DICRT(EL=dat$E_L, ER=dat$E_R, SL=dat$S_L, SR=dat$S_R, iter=1000)
plot(out$value, out$prob, type="s", ylab="eccdf", xlab="x")
curve(pweibull(x,1.5,5,lower.tail=FALSE), add=TRUE, col="royalblue")
}