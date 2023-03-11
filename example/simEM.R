####
set.seed(999)
n <- 100L
x <- rweibull(n,1.5,5)
dat <- simDIC(x)
out_em <- eccdf_dic_em(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                       ctype=rep(3L,n), iter=500L)
sum(out_em$event)
plot(out_em$lp, type = "l")

ind=with(out_em, event>0)
vv <-with(out_em, (prob^2)*(event-prob^2*sum(event)*(1-prob)))

z <- qnorm(0.975)
b <- with(out_em, rev(exp(z*sqrt(cumsum(rev(vv)))/log(ccdf))))
dfout <- data.frame(value=out_em$value,
                    ccdf=out_em$ccdf[-length(out_em$ccdf)],
                    b=b[-length(out_em$ccdf)])

ggplot(dfout, aes(x=value, y=ccdf)) +
  geom_line()+
  geom_ribbon(aes(ymax=ccdf^b, ymin=ccdf^(1/b)), alpha=0.1)+
  stat_function(fun = pweibull, args = list(shape=1.5, scale=5, lower.tail=FALSE), colour="royalblue")+
  labs(y="ccdf", x="y")+
  theme_classic(16)

