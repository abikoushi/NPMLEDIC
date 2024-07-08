library(NPMLEDIC)
set.seed(222); y <- rweibull(50,1.5,10)
dat <- simDIC(y)
breaks <- sort(unique(c(dat$LE, dat$RE, dat$LS, dat$RS)))
LE <- match(dat$LE,breaks)
RE <- match(dat$RE,breaks)
LS <- match(dat$LS,breaks)
RS <- match(dat$RS,breaks)

res <- NPMLEDIC:::joint_DIC_em(LE, RE, LS, RS, 
                               m = length(breaks),
                               maxit=2, tol=0.01)
#image(res$event)
plot(res$lp)
plot(breaks, cumsum(res$alpha)/sum(res$alpha), type = "s")
curve(pweibull(x,1.5,10),add=TRUE, col="royalblue")

plot(breaks, cumsum(res$beta), type = "s")
abline(0,1,col="royalblue")
