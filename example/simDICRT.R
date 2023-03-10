library(NPMLEDIC)
library(ggplot2)

set.seed(505)
n <- 25L
x <- rweibull(n,1.5,5)
dat <- simDIC(x)
out <- eccdf_dic_vb(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                    ctype=rep(3L,n), alpha0 = 0.01, iter=500L)
sum(out$event)
plot(out$lp, type = "l")

ss <- reccdf(2000, out$alpha)
ci <- apply(ss, 1, quantile, prob=c(0.975,0.025))

# with(out, plot(value, ccdf[-length(ccdf)], type="s"))
# with(out, lines(value, ci[1,-length(ccdf)], type="s", lty=2))
# with(out, lines(value, ci[2,-length(ccdf)], type="s", lty=2))
# curve(pweibull(x, shape=1.5, scale=5, lower.tail=FALSE), add=TRUE, col="royalblue")

dfout <- data.frame(value=out$value,
                    ccdf=out$ccdf[-length(out$ccdf)],
                    upper=ci[1,-length(out$ccdf)],
                    lower=ci[2,-length(out$ccdf)])
ggplot(dfout, aes(x=value, y=ccdf)) +
  geom_step()+
  geom_ribbon(aes(ymax=upper, ymin=lower), alpha=0.1)+
  stat_function(fun = pweibull, args = list(shape=1.5, scale=5, lower.tail=FALSE), colour="royalblue")+
  labs(y="ccdf", x="y")+
  theme_classic(16)
ggsave("CI.png")

# ind=with(out, event>0)
# vv <-with(out, (prob[ind]^2)*(event[ind]-prob[ind]^2*sum(event[ind])*(1-prob[ind])))
# with(out, plot(value[ind], ccdf[ind], type="s"))
# z <- qnorm(0.975)
# b <- with(out, rev(exp(z*sqrt(cumsum(rev(vv)))/log(ccdf[ind]))))

dfout <- data.frame(value=out$value[ind], ccdf=out$ccdf[ind])
ggplot(dfout, aes(x=value, y=ccdf)) +
  geom_step()+
  geom_ribbon(aes(ymax=ccdf^b, ymin=ccdf^(1/b)), alpha=0.1)+
  stat_function(fun = pweibull, args = list(shape=1.5, scale=5, lower.tail=FALSE), colour="royalblue")+
  labs(y="ccdf", x="y")+
  theme_classic(16)

