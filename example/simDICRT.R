library(NPMLEDIC)
library(ggplot2)

set.seed(457)
n <- 25L
x <- rweibull(n,1.5,5)
dat <- simDIC(x)
out <- eccdf_dic_vb(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                    ctype=rep(3L,n), alpha0=0.5, iter=100L)

plot(out$lp,type = "l")

with(out, plot(value, ccdf, type="s"))
curve(pweibull(x, shape=1.5, scale=5, lower.tail=FALSE), add=TRUE, col="royalblue")

randp <- reccdf(2000, out$alpha)
ci <- apply(randp, 1, quantile, prob=c(0.025,0.975))
p <- pweibull(out$value, shape=1.5, scale=5, lower.tail=FALSE)
mean(ci[1,]<=p & p<=ci[2,])

ggplot(out, aes(x=value, y=prob)) +
  geom_step()+
  geom_ribbon(aes(ymin=ci[1,], ymax=ci[2,]), alpha=0.1)+
  stat_function(fun = pweibull, args = list(shape=1.5, scale=5, lower.tail=FALSE), colour="royalblue")+
  labs(y="eccdf", x="y")+
  theme_classic(16)

