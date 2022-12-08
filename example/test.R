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

#####
res <- vector("list", 100)
for(i in 1:100){
  x <- rweibull(100,2,7)
  dat <- simDIC(x)
  out <- eccdf_DIC(dat$LE, dat$RE, dat$S, dat$S, rep(1L,100), iter = 10)
  res[[i]] <- data.frame(out,id=i)
}
# with(out,plot(value, prob, type="s"))
# curve(pweibull(x,2,7,lower.tail = FALSE),add=TRUE, lty=2)
dfs <- bind_rows(res)

ggplot(dfs, aes(x=value, y=prob))+
  geom_step(aes(group=id), alpha=0.05)+
  stat_function(fun=function(x)pweibull(x,2,7, lower.tail = FALSE), colour="royalblue")+
  theme_bw()

ggsave("test.png")
