set.seed(123456)
n <- 100L
x <- rweibull(n,1.5,5)
dat <- simDIC(x)
out <- eccdf_DIC(EL=dat$LE, ER=dat$RE, SL=dat$LS, SR=dat$RS,
                 ctype=rep(3L,n), iter=1000)
with(out, plot(value, prob, type="s", ylab="eccdf", xlab="x"))
curve(pweibull(x,1.5,5,lower.tail=FALSE), add=TRUE, col="royalblue")

set.seed(12345)
n <- 100L
at <- runif(n, 0, 4)
x <- rweibull(n,1.5,5)
dat0 <- simDIC(x, at = at)
dat <- dat0[dat0$RS<=7,]
dim(dat)
out <- eccdf_DICT(EL=dat$LE, ER=dat$RE, SL=dat$LS, SR=dat$RS,
                  tmax = rep(6,nrow(dat)),
                  ctype=rep(3L,nrow(dat)), iter=1000)
with(out, plot(value, prob, type="s", ylab="eccdf", xlab="x"))
curve(pweibull(x,1.5,5,lower.tail=FALSE), add=TRUE, col="royalblue")

