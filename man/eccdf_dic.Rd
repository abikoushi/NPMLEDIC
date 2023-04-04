\name{eccdf_dic}
\alias{eccdf_dic_em}
\alias{eccdf_dic_vb}
\alias{eccdf_dic_gibbs}

\title{
Estimate emperical survival function from doubly interval censored data
}
\description{
Estimate emperical survival function from doubly interval censored  data
}
\usage{
eccdf_dic_em(EL, ER, SL, SR, ctype, iter=1000L)

eccdf_dic_vb(EL, ER, SL, SR, ctype,
          alpha0 = 1, iter=1000L)
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
\item{alpha0}{
parameter of the dirichlet prior
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
out <- eccdf_dic_vb(EL=dat$LE, ER=dat$RE, SL=dat$LS, SR=dat$RS,ctype=rep(3L,n), iter=1000)
with(out, plot(value, rmlast(ccdf), type="s", ylab="eccdf", xlab="x"))
curve(pweibull(x,1.5,5,lower.tail=FALSE), add=TRUE, col="royalblue")
}