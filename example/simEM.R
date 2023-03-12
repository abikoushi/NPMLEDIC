library(NPMLEDIC)
library(ggplot2)

set.seed(888)
n <- 50L
x <- rweibull(n,1.5,5)
dat <- simDIC(x)
out_em <- eccdf_dic_em(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                       ctype=rep(3L,n), alpha0 = 0.1, iter=500L)
sum(out_em$event)
plot(out_em$lp, type = "l")
names(out_em)

dfout <- confint_dic(out_em)

ggplot(dfout, aes(x=value, y=ccdf)) +
  geom_step()+
  geom_step(aes(y=lower),alpha=0.3)+
  geom_step(aes(y=upper),alpha=0.3)+
  stat_function(fun = pweibull, args = list(shape=1.5, scale=5, lower.tail=FALSE), colour="royalblue")+
  labs(y="ccdf", x="y")+
  theme_classic(16)

