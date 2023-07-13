Example: Flu A
================
Ko Abe
2023-07-13

``` r
library(NPMLEDIC)
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(coarseDataTools)
library(pammtools)
```

    ## 
    ## Attaching package: 'pammtools'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
data(fluA.inc.per)
head(fluA.inc.per)
```

    ##   author year EL  ER  SL  SR type
    ## 1  moser 1977  0 0.5 1.0 1.5    0
    ## 2  moser 1977  0 0.5 1.5 2.0    0
    ## 3  moser 1977  0 0.5 1.5 2.0    0
    ## 4  moser 1977  0 0.5 1.5 2.0    0
    ## 5  moser 1977  0 0.5 1.5 2.0    0
    ## 6  moser 1977  0 0.5 1.5 2.0    0

``` r
fluAdf <- mutate(fluA.inc.per,ctype=type) %>% 
  mutate(ctype=if_else(SL==SR & EL==ER, 0L, 3L)) %>% 
  mutate(ctype=if_else(SL==SR & EL<ER, 1L, ctype)) %>% 
  mutate(ctype=if_else(SL<SR & EL==ER, 2L, ctype)) 
```

``` r
out_flu <- eccdf_dic_em(LE = fluAdf$EL, RE = fluAdf$ER,
                        LS = fluAdf$SL, RS = fluAdf$SR,
                        ctype = fluAdf$ctype, alpha0 = 1e-8)
```

``` r
probs=c(0.95)
ci_flu <- confint_dic(out_flu, prob=probs)
```

``` r
flu_opt_w <- optim(c(0,1),loglikhd,dat=fluA.inc.per,dist = "W")
flu_opt_w <- optim(flu_opt_w$par,loglikhd,dat=fluA.inc.per,dist = "W",
                   method = "BFGS")
flu_opt_g <- optim(c(0,1),loglikhd,dat=fluA.inc.per,dist = "G")
flu_opt_g <- optim(flu_opt_g$par,loglikhd,dat=fluA.inc.per,dist = "G",
                   method = "BFGS")
flu_opt_l <- optim(c(0,1),loglikhd,dat=fluA.inc.per,dist = "L")
flu_opt_l <- optim(flu_opt_l$par,loglikhd,dat=fluA.inc.per,dist = "L",
                   method = "BFGS")
```

``` r
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
```

![](fluA_files/figure-gfm/plot%20results-1.png)<!-- -->
