library(NPMLEDIC)
library(ggplot2)
library(dplyr)
library(coarseDataTools)
library(pammtools)

#Flu A
data(fluA.inc.per)
head(fluA.inc.per)
fluAdf <- mutate(fluA.inc.per,ctype=type) %>% 
  mutate(ctype=if_else(SL==SR & EL==ER, 0L, 3L)) %>% 
  mutate(ctype=if_else(SL==SR & EL<ER, 1L, ctype)) %>% 
  mutate(ctype=if_else(SL<SR & EL==ER, 2L, ctype)) 

we <- dplyr::filter(fluAdf, ctype!=0) %>% 
  mutate(WE=ER-EL,type="E")
ws <- dplyr::filter(fluAdf, ctype==3) %>% 
  mutate(WS=SR-SL,type="S")

tab <- table(fluAdf$ctype)
xtable(t(tab))

df <- data.frame(w=c(we$WE,ws$WS), censor=c(we$type,ws$type))
sd(we$WE)
sd(ws$WS)
ggplot(df, aes(x=w, linetype=censor))+
  stat_ecdf()+
  labs(x="width")+
  theme_bw(16)

out_flu <- eccdf_dic_em(LE = fluAdf$EL, RE = fluAdf$ER,
                        LS = fluAdf$SL, RS = fluAdf$SR,
                        ctype = fluAdf$ctype, alpha0 = 1e-8)


probs=c(0.95)
ci_flu <- confint_dic(out_flu, prob=probs)

flu_opt_w <- optim(c(0,1),loglikhd,dat=fluA.inc.per,dist = "W")
flu_opt_w <- optim(flu_opt_w$par,loglikhd,dat=fluA.inc.per,dist = "W",
                   method = "BFGS")
flu_opt_g <- optim(c(0,1),loglikhd,dat=fluA.inc.per,dist = "G")
flu_opt_g <- optim(flu_opt_g$par,loglikhd,dat=fluA.inc.per,dist = "G",
                   method = "BFGS")
flu_opt_l <- optim(c(0,1),loglikhd,dat=fluA.inc.per,dist = "L")
flu_opt_l <- optim(flu_opt_l$par,loglikhd,dat=fluA.inc.per,dist = "L",
                   method = "BFGS")


ggplot(ci_flu, aes(x=value))+
  geom_stepribbon(aes(ymin=lower, ymax=upper), alpha=0.3)+
  geom_step(aes(y=ccdf))+
  stat_function(fun=pgamma, args = list(shape=exp(flu_opt_g$par[1]), scale=exp(flu_opt_g$par[2]), lower.tail=FALSE),
                mapping = aes(linetype="gamma", colour="gamma"))+
  stat_function(fun=pgamma, args = list(shape=exp(flu_opt_g$par[1]), scale=exp(flu_opt_g$par[2]), lower.tail=FALSE),
                geom = "point", n=20, aes(shape="gamma", colour="gamma"), size=2.5)+
  stat_function(fun=plnorm, args = list(meanlog=flu_opt_l$par[1], sdlog=exp(flu_opt_l$par[2]), lower.tail=FALSE),
                aes(linetype="log normal", colour="log normal"))+
  stat_function(fun=plnorm, args = list(meanlog=flu_opt_l$par[1], sdlog=exp(flu_opt_l$par[2]), lower.tail=FALSE),
                geom = "point", n=20, aes(shape="log normal", colour="log normal"), size=2.5)+
  stat_function(fun=pweibull, args = list(shape=exp(flu_opt_w$par[1]), scale=exp(flu_opt_w$par[2]), lower.tail=FALSE),
                aes(linetype="weibull", colour="weibull"))+
  stat_function(fun=pweibull, args = list(shape=exp(flu_opt_w$par[1]), scale=exp(flu_opt_w$par[2]), lower.tail=FALSE),
                geom = "point", n=20, aes(shape="weibull", colour="weibull"), size=2.5)+
  labs(x="incubation period", y="ccdf", shape="model", linetype="model", colour="model")+
  scale_colour_manual(values = c("orange2", "royalblue", "forestgreen"))+
  theme_bw(16)
