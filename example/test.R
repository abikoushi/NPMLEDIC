library(NPMLEDIC)
library(ggplot2)
library(dplyr)
library(readr)
library(ggplot2)
library(pammtools)

integrate(function(x)-dunif(x,-1.5,1.5, log = TRUE)*dnorm(x), -1.5, 1.5)
integrate(function(x)-dunif(x,-0.5,0.5, log = TRUE)*dnorm(x), -0.5, 0.5)

dat <- read_csv("./example/dat.csv")
dat_t <- dplyr::select(dat, EL,ER,SL,SR) %>% 
  dplyr::filter(!is.na(EL)) %>% 
  mutate(startrec = min(EL,ER,SL,SR)) %>% 
  mutate(EL=as.integer(EL-startrec),ER=as.integer(ER-startrec),
         SL=as.integer(SL-startrec),SR=as.integer(SR-startrec)) %>% 
  mutate(ctype=if_else(SL==SR & EL==ER, 0L, 1L)) %>% 
  mutate(ctype=if_else(SL==SR & EL<ER, 1L, ctype)) %>% 
  mutate(ctype=if_else(SL<SR & EL==ER, 2L, ctype)) %>% 
  mutate(ctype=if_else(SL<SR & EL<ER, 3L, ctype)) 

out <- eccdf_dic_em(dat_t$EL, dat_t$ER, dat_t$SL, dat_t$SR, dat_t$ctype,
                    alpha0 = 1e-9)

plot(out$lp,type = "l")

probs=c(0.6,0.7,0.8,0.95)
out_ci <- confint_dic(out, prob=probs, scale = "linear")

ggplot(out_ci$confint, aes(x=value))+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=level, group=level), alpha=0.2)+
  geom_line(data=out_ci$point,aes(y=ccdf))+
  labs(fill="level")+
  theme_bw(16)
#ggsave("ci.png")
##
library(survival)
w_shape <- 1.5
w_scale <- 7

g_shape <- 1.5
g_scale <- 4
#
cat("mean: ",gamma(1+1/w_shape)*w_scale)
cat("mean: ",g_shape*g_scale)
#sfit <- survfit(Surv(time = dat$S-dat$RE, time2 = dat$S-dat$LE, type="interval2")~1)
# plot(sfit, col="royalblue")
# with(df,lines(x = value, y=ccdf), type="s")
# with(df,lines(x = value, y=upper, type="s", lty=2))
# with(df,lines(x = value, y=lower, type="s", lty=2))


x <- rweibull(500,w_shape,w_scale)
#x <- rgamma(100,g_shape,scale=g_scale)
dat <- simDIC(x)

out <- eccdf_dic_em(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                       alpha0 = 1e-10, ctype = 3L, maxit = 5000)
length(out$value)
length(out$prob)
plot(out$lp,type="l")
probs=c(0.8, 0.9, 0.95, 0.99)
out_ci <- confint_dic(out, prob = probs, scale = "loglog")
out_ci <- confint_dic(out, prob = probs, scale = "linear")
ggplot(out_ci$confint, aes(x=value))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=level, group = reorder(level, -level)), alpha=0.2)+
  geom_line(data=out_ci$point, aes(y=ccdf))+
  #stat_function(fun=pgamma, args = list(shape=g_shape,scale=g_scale,lower.tail=FALSE), colour="firebrick")+
  stat_function(fun=pweibull, args = list(shape=w_shape,scale=w_scale,lower.tail=FALSE), colour="firebrick")+
  scale_fill_gradient(low="gray10",high = "gray90")+
  theme_bw(16)

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

