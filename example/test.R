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

####matrixでalpha持つ
library(ggplot2)
x <- rweibull(100,2,7.5)
dat <- simDIC(x)
dat <- dat[dat$RS<100,]
nrow(dat)
out <- eccdf_dicrt_em(LE=dat$LE, RE=dat$RE, LS=dat$LS, RS=dat$RS,
                      alpha0 = 1e-8,
                      tmax = 100, ctype = 3L, iter = 500)
plot(out$lp,type="l")
sum(out$event)
sum(out$event2)
df <- confint_dic(out)
ggplot(df, aes(x=value, y=ccdf))+
  geom_step()+
  geom_step(aes(y=lower), alpha=0.4)+
  geom_step(aes(y=upper), alpha=0.4)+
  stat_function(fun=pweibull, args = list(shape=2,scale=7.5,lower.tail=FALSE), colour="royalblue")+
  theme_bw(16)
