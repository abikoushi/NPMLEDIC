\name{eccdf_DIC}
\alias{eccdf_DIC}
\title{
Estimate emperical survival function from doubly interval censored data
}
\description{
Estimate emperical survival function from doubly interval censored  data
}
\usage{
eccdf_DIC(EL, ER, SL, SR, breaks=NULL, iter=1000L)
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
Estimate emperical survival function from doubly interval censored data. 
}
\examples{
set.seed(123456)
n <- 100L
x <- rweibull(n,1.5,5)
dat <- simDIC(x, ew, sw)
out <- eccdf_DIC(EL=dat$E_L, ER=dat$E_R, SL=dat$S_L, SR=dat$S_R, iter=1000)
with(out, plot(value, prob, type="s", ylab="eccdf", xlab="x"))
curve(pweibull(x,1.5,5,lower.tail=FALSE), add=TRUE, col="royalblue")
}