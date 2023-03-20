library(NPMLEDIC)
library(ggplot2)
library(dplyr)
library(readr)
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

out <- eccdf_DIC(dat_t$EL, dat_t$ER, dat_t$SL, dat_t$SR, dat_t$ctype, iter = 100)
with(out,plot(value,prob,type="s"))

library(ggplot2)
library(pammtools)
w_shape <- 2
w_scale <- 5
gamma(1+1/w_shape)*w_scale
x <- rweibull(100,w_shape,w_scale)
dat <- simDIC(x)
out <- eccdf_dic_em(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                       alpha0 = 0.1, ctype = 3L, iter = 250)

plot(out$lp,type="l")

sd <- rev(sqrt(cumsum(rev(diag(solve(out$I))))))
sd <- sd[-length(sd)]
df <- confint_dic(out,0.95)

# df <- data.frame(value=out$value, ccdf=rmlast(out$ccdf),
#                  sd=sd)
ggplot(df, aes(x=value, y=ccdf))+
  geom_step()+
  geom_stepribbon(aes(ymin=lower,ymax=upper), alpha=0.3)+
  stat_function(fun=pweibull, args = list(shape=w_shape,scale=w_scale,lower.tail=FALSE), colour="royalblue")+
  theme_bw(16)

########

out <- eccdf_dic_gibbs(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                    alpha0 = 1, ctype = 3L, iter = 2500)
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



# ggplot(df, aes(x=value, y=ccdf))+
#   geom_step()+
#   geom_step(aes(y=lower), alpha=0.4)+
#   geom_step(aes(y=upper), alpha=0.4)+
#   stat_function(fun=pweibull, args = list(shape=2,scale=7.5,lower.tail=FALSE), colour="royalblue")+
#   theme_bw(16)



#######
dat <- dat[dat$RS<100,]
nrow(dat)
out <- eccdf_dicrt_em(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                      alpha0 = 1e-8,
                      tmax = 100, ctype = 3L, iter = 500)
plot(out$lp,type="l")
sum(out$event)
sum(out$event2)
