\name{eccdf_DIC}
\alias{eccdf_DIC}
\alias{eccdf_DICT}

\title{
Estimate emperical survival function from doubly interval censored data
}
\description{
Estimate emperical survival function from doubly interval censored  data
}
\usage{
eccdf_DIC(EL, ER, SL, SR, ctype,
          breaks=NULL, iter=1000L)

eccdf_DICT(EL, ER, SL, SR, ctype,
           tmin=-Inf, tmax=Inf,
           breaks=NULL, iter=1000L)
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
breaks of histrogram-type estimator
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
dat <- simDIC(x)
out <- eccdf_DIC(EL=dat$LE, ER=dat$RE, SL=dat$LS, SR=dat$RS,
                 ctype=rep(3L,n), iter=1000)
with(out, plot(value, prob, type="s", ylab="eccdf", xlab="x"))
curve(pweibull(x,1.5,5,lower.tail=FALSE), add=TRUE, col="royalblue")
}