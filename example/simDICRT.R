library(NPMLEDIC)
library(ggplot2)

set.seed(123456)
n <- 100L
x <- rweibull(n,1.5,5)
dat <- simDIC(x)
out <- eccdf_DIC(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                 ctype=rep(3L,n), iter=1000)
with(out, plot(value, prob, type="s", ylab="eccdf", xlab="x"))
curve(pweibull(x,1.5,5,lower.tail=FALSE), add=TRUE, col="royalblue")

set.seed(12345)
n <- 100L
at <- runif(n, 0, 4)
x <- rweibull(n,1.5,5)
dat0 <- simDIC(x, at = at)
dat0$trunc <- dat0$RS>7

ggplot(dat0, aes(x=LE,fill=trunc))+
  geom_histogram()

ggplot(dat0, aes(x=RE,fill=trunc))+
  geom_histogram()

out <- eccdf_DICT(EL=dat$LE, ER=dat$RE, SL=dat$LS, SR=dat$RS,
                  tmax = rep(6,nrow(dat)),
                  ctype=rep(3L,nrow(dat)), iter=1000)
with(out, plot(value, prob, type="s", ylab="eccdf", xlab="x"))
curve(pweibull(x,1.5,5,lower.tail=FALSE), add=TRUE, col="royalblue")

