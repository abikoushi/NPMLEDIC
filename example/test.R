library(NPMLEDIC)
library(ggplot2)
library(dplyr)
library(readr)
library(pammtools)
library(survival)

g_shape <- 2
g_scale <- 5
cat("Gamma ","mean: ",g_shape*g_scale, " variance: ", g_shape*g_scale^2)

w_shape <- 1.5
w_scale <- 10/gamma(1+1/w_shape)
#w_scale <- 5/sqrt((gamma(1+2/w_shape) - gamma(1+1/w_shape)^2))
cat("Weibull ", "mean: ",gamma(1+1/w_shape)*w_scale,  " variance: ", (w_scale^2)*(gamma(1+2/w_shape) - gamma(1+1/w_shape)^2))

sigma <- 0.5
mu <- log(10/exp(sigma^2/2))
cat("LogNormal ", "mean: ",exp(mu+sigma^2/2), " variance: ", (exp(sigma^2)-1)*exp(2*mu+sigma^2))

n <- 50L
#x <- rweibull(n, w_shape, w_scale)
x <- rgamma(n, g_shape, scale=g_scale)
#x <- rlnorm(n, mu, sigma)
dat <- simDIC(x, WS = 2, WE = 2)
system.time({
  out_em <- eccdf_dic_em(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                      ctype = 3L, maxit = 1000L)
})

#
system.time({
  out_vb <- eccdf_dic_vb(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                      ctype = 3L, maxit  = 1000L)
})

plot(out_vb$lp, type="l")
probs=c(0.95)
head(out$estimates)
#out_ci <- confint_dic(out, prob = probs, scale = "loglog")
out_ci <- confint_dic(out_em, prob = probs)

ggplot(out_ci, aes(x=value))+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=level, group = reorder(level, -level)), alpha=0.2)+
  geom_line(data=out_em$estimates, aes(y=ccdf))+
  stat_function(fun=pgamma, args = list(shape=g_shape,scale=g_scale,lower.tail=FALSE), colour="orange2")+
  #stat_function(fun=pweibull, args = list(shape=w_shape,scale=w_scale,lower.tail=FALSE), colour="orange2")+
  #stat_function(fun=plnorm, args = list(meanlog=mu, sdlog=sigma,lower.tail=FALSE), colour="orange2")+
  #scale_fill_gradient(low="gray10", high = "gray90")+
  theme_bw(16)
head(out_ci)
ggplot(out_ci, aes(x=pgamma(value,shape=g_shape,scale=g_scale,lower.tail=FALSE)-ccdf, y=se, colour=value))+
  geom_point()+
  facet_wrap(~level)+
  theme_bw(16)


rand <- reccdf(5000L, out_vb$estimates$alpha)
ci <- apply(rand, 1, quantile, prob=c(0.025,0.975))
pmean <- apply(rand, 1, mean)

plot(out_vb$estimates$value,pmean,type="s")
curve(pgamma(x, g_shape, scale=g_scale, lower.tail = FALSE), add=TRUE, col="royalblue")
lines(out_vb$estimates$value, ci[1,], col="royalblue", lty=2)
lines(out_vb$estimates$value, ci[2,], col="royalblue", lty=2)
