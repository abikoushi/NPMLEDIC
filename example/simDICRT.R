library(NPMLEDIC)
library(ggplot2)
library(pammtools)
library(parallel)

qmix <- function(p,w,a,b,upper=100){
  cdf1 <- function(x)pweibull(x,a,b)
  cdf2 <- function(x)pexp(x)
  sol <- sapply(p, function(p)uniroot(function(x){w*cdf1(x)+(1-w)*cdf2(x)-p}, c(0, upper=upper))$root)
  return(sol)
}

rmix <- function(n,w,a,b){
  m <- rbinom(1,n,w)
  x1 <- rweibull(m,a,b)
  x2 <- rexp(n-m)
  sample(c(x1,x2))
}

pmix <- function(x, w, a, b){
  ccdf1 <- function(x)pweibull(x, a, b, lower.tail = FALSE)
  ccdf2 <- function(x)pexp(x, lower.tail = FALSE)
  w*ccdf1(x)+(1-w)*ccdf2(x)
}

# curve(pmix(x,0.8,4,10),0,25)
# x <- rmix(5000,0.8,4,10)
# mean(x)
# hist(x, breaks = "FD")
# q <- qmix(0.5,0.8,4,10)
# # mean(x<q) # test median

qeccdf <- function(p0, p, breaks){
  sapply(p0, function(x){f=1-p<=x;
  ifelse(all(f),Inf,breaks[which.min(f)])})
}

simcp <- function(i,n){
  x <- rmix(n,0.8,4,10)
  dat <- simDIC(x)
  out_vb <- eccdf_dic_vb(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                         ctype=rep(3L,n), alpha0 = 1e-8, iter=500L)
  ss <- reccdf(2000, out_vb$alpha)
  ci_lower <- apply(ss, 1, quantile, prob=c(0.1,0.05,0.025,0.005))
  ci_upper <- apply(ss, 1, quantile, prob=c(0.9,0.95,0.975,0.995))
  q_true <- qmix(c(0.5,0.8,0.9,0.95,0.99), w=0.8,a=4,b=10)
  q_low <- apply(ci_lower, 1, function(x)qeccdf(c(0.5,0.8,0.9,0.95,0.99), x, out_vb$value))
  q_up <- apply(ci_upper, 1, function(x)qeccdf(c(0.5,0.8,0.9,0.95,0.99), x, out_vb$value))
  cp <- q_low <= q_true & q_true <= q_up
  colnames(cp) <- c("0.8", "0.9", "0.95", "0.99")
  row.names(cp) <- c("0.5", "0.8", "0.9", "0.95", "0.99")
  dfout <- data.frame(value=out_vb$value,
                      ccdf=rmlast(out_vb$ccdf))
  list(ccdf=dfout, cp=cp)
}
set.seed(999)
simres <- mclapply(1:100, simcp, n=100, mc.cores=detectCores())
apply(simplify2array(lapply(simres, function(x)x$cp)), 1:2, mean)

dfsimres <- do.call("rbind", lapply(1:100, function(i)cbind(simres[[i]]$ccdf,id=i)))

ggplot(dfsimres, aes(x=value, y=ccdf)) +
  geom_step(aes(group=id),alpha=0.05)+
  stat_function(fun = pmix, args = list(w=0.8,a=4,b=10),colour="blue")+
  labs(y="ccdf", x="y")+
  theme_classic(16)

