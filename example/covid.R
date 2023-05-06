library(NPMLEDIC)
library(ggplot2)
library(dplyr)
library(readr)
library(pammtools)
library(survival)

# integrate(function(x)-dunif(x,-1.5,1.5, log = TRUE)*dnorm(x), -1.5, 1.5)
# integrate(function(x)-dunif(x,-0.5,0.5, log = TRUE)*dnorm(x), -0.5, 0.5)

# dat <- read_csv("./example/dat.csv")
# dat_t <- dplyr::select(dat, EL,ER,SL,SR) %>% 
#   dplyr::filter(!is.na(EL)) %>% 
#   mutate(startrec = min(EL,ER,SL,SR)) %>% 
#   mutate(EL=as.integer(EL-startrec),ER=as.integer(ER-startrec),
#          SL=as.integer(SL-startrec),SR=as.integer(SR-startrec)) %>% 
#   mutate(ctype=if_else(SL==SR & EL==ER, 0L, 1L)) %>% 
#   mutate(ctype=if_else(SL==SR & EL<ER, 1L, ctype)) %>% 
#   mutate(ctype=if_else(SL<SR & EL==ER, 2L, ctype)) %>% 
#   mutate(ctype=if_else(SL<SR & EL<ER, 3L, ctype)) 

# out <- eccdf_dic_em(dat_t$EL, dat_t$ER, dat_t$SL, dat_t$SR, dat_t$ctype,
#                     alpha0 = 1e-9)

# plot(out$lp,type = "l")
# 
# probs=c(0.6,0.7,0.8,0.95)
# out_ci <- confint_dic(out, prob=probs, scale = "linear")
# 
# ggplot(out_ci$confint, aes(x=value))+
#   geom_ribbon(aes(ymin=lower, ymax=upper, fill=level, group=level), alpha=0.2)+
#   geom_line(data=out_ci$point,aes(y=ccdf))+
#   labs(fill="level")+
#   theme_bw(16)
#ggsave("ci.png")
##