library(NPMLEDIC)
library(ggplot2)
library(dplyr)
library(readr)
library(pammtools)
library(survival)

w_shape <- 1.5
w_scale <- 7

g_shape <- 1.5
g_scale <- 4

mu <- 1.7
sigma <- 0.5

cat("mean: ",exp(mu+sigma^2/2), " variance: ", (exp(sigma^2)-1)*exp(2*mu+sigma^2))
cat("mean: ",gamma(1+1/w_shape)*w_scale,  " variance: ", (w_scale^2)*(gamma(1+2/w_shape) - gamma(1+1/w_shape)^2))
cat("mean: ",g_shape*g_scale, " variance: ", g_shape*g_scale^2)

n <- 500L
x <- rweibull(n, w_shape, w_scale)
#x <- rgamma(n, g_shape, scale=g_scale)
#x <- rlnorm(n, mu, sigma)
dat <- simDIC(x)
system.time({
  out <- eccdf_dic_em(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                      alpha0 = 1e-5, ctype = 3L, maxit = 10000L)
})
head(out$value)
plot(out$lp,type="l")
probs=c(0.8, 0.9, 0.95, 0.99)
#out_ci <- confint_dic(out, prob = probs, scale = "loglog")
out_ci <- confint_dic(out, prob = probs, scale = "linear")
ggplot(out_ci$confint, aes(x=value))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=level, group = reorder(level, -level)), alpha=0.2)+
  geom_line(data=out_ci$point, aes(y=ccdf))+
  #stat_function(fun=pgamma, args = list(shape=g_shape,scale=g_scale,lower.tail=FALSE), colour="orange2")+
  stat_function(fun=pweibull, args = list(shape=w_shape,scale=w_scale,lower.tail=FALSE), colour="orange2")+
  #stat_function(fun=plnorm, args = list(meanlog=mu, sdlog=sigma,lower.tail=FALSE), colour="orange2")+
  #scale_fill_gradient(low="gray10", high = "gray90")+
  theme_bw(16)

###
n <- 100L
#x <- rweibull(n, w_shape, w_scale)
x <- rgamma(n, g_shape, scale=g_scale)
dat <- simDIC(x)
breaks <- with(dat, sort(unique(c(0, S-RE, S-LE))))
breaks <- breaks[breaks>=0]
breaks
system.time({
  out <- eccdf_dic_em(LE=dat$LE, RE=dat$RE, LS=dat$S, RS=dat$S,
                      breaks = breaks,
                      alpha0 = 1e-8, ctype = 1L, maxit = 300, tol = 1e-8)
})
ind=with(out, event>0)
variance <- diag(solve(out$I[ind,ind]))
z <- qnorm(1-0.5*(1-probs))
value = out$value[ind]
ccdf = out$ccdf[ind]
if(scale=="linear"){
  b <- sapply(z, function(z0)rev(z0*sqrt(cumsum(rev(variance)))))
  lower = as.matrix(ifelse(ccdf-b<0, 0, ccdf-b))
  upper = as.matrix(ifelse(ccdf+b>1, 1, ccdf+b))
}

matplot(out$value, upper, type="l", lty=1, col=grey.colors(4, start=0.5, end=0.8))
matplot(out$value, lower, type="l", lty=1, col=grey.colors(4, start=0.5, end=0.8), add=TRUE)
with(out, lines(value, ccdf[-length(ccdf)], type="s"))
curve(pgamma(x, shape = 1.5, scale = 4, lower.tail = FALSE), add=TRUE, col="royalblue")
sum(out$event[-length(out$event)])
#out$event[length(out$event)]
length(out$prob)
#out$lp
plot(out$lp,type="l")




########

out <- eccdf_dic_gibbs(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                    alpha0 = 0.1, ctype = 3L, iter = 2500)
plot(out$lp,type="l")
burnin <- 1:500
mccdf <- rowMeans(out$ccdf[-nrow(out$ccdf),-burnin])
qccdf <- apply(out$ccdf[-nrow(out$ccdf),-burnin],1,quantile, prob=c(0.025, 0.975))
dim(out$ccdf[-nrow(out$ccdf),-burnin])
df <- data.frame(value = out$value,
                 ccdf=mccdf,
                 lower=qccdf[1,],
                 upper=qccdf[2,])

ggplot(df, aes(x=value, y=ccdf))+
  geom_step()+
  geom_stepribbon(aes(ymin=lower, ymax=upper),alpha=0.2)+
  stat_function(fun=pweibull, args = list(shape=1.5,scale=8,lower.tail=FALSE), colour="royalblue")+
  theme_bw(16)

