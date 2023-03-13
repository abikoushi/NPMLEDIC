
p2ccdf <- function(p){
  rev(cumsum(rev(p)))
}

#' @export rmlast
rmlast <- function(x)x[-length(x)]

#' @export eccdf_dic_em
eccdf_dic_em <- function(LE, RE, LS, RS, ctype, alpha0 = 0, iter=1000L){
  breaks <- sort(unique(c(0,RS-RE, RS-LE, LS-RE, LS-LE)))
  breaks <- breaks[breaks>=0]
  
  n <- length(LE)
  if(length(ctype)==1){
    ctype = rep(ctype, n)
  }
  res <- ep_DIC_em(LE, RE, LS, RS, ctype, breaks, alpha0, iter)    
  res$value <- breaks
  res$ccdf  <- with(res, p2ccdf(prob))
  res$variance <-with(res, (prob^2)/(event-prob^2*sum(event)*(1-prob)))
  return(res)
}

#' @export eccdf_dic_em
eccdf_dicrt_em <- function(LE, RE, LS, RS, tmax, ctype, alpha0 = 0, iter=1000L){
  breaks <- sort(unique(c(0,RS-RE, RS-LE, LS-RE, LS-LE)))
  breaks <- breaks[breaks>=0]
  n <- length(LE)
  if(length(ctype)==1L){
    ctype = rep(ctype, n)
  }
  if(length(tmax)==1L){
    tmax = rep(tmax, n)
  }
  res <- ep_DICT_em(LE, RE, LS, RS, tmax, ctype, breaks, alpha0, iter)    
  res$value <- breaks
  res$ccdf  <- with(res, p2ccdf(prob))
  res$variance <-with(res, (prob^2)/(event-prob^2*sum(event)*(1-prob)))
  return(res)
}

#' @export eccdf_dic_vb
eccdf_dic_vb <- function(LE, RE, LS, RS, ctype, alpha0 = 1, iter = 1000L){
  breaks <- sort(unique(c(0,RS-RE, RS-LE, LS-RE, LS-LE)))
  breaks <- breaks[breaks>=0]
  
  n <- length(LE)
  if(length(ctype)==1L){
    ctype = rep(ctype, n)
  }
  
  res <- ep_DIC_vb(LE, RE, LS, RS, ctype, breaks, alpha0, iter)
  res$value <- breaks
  res$ccdf  <- with(res, p2ccdf(alpha/sum(alpha)))
  return(res)
}

#' @export reccdf
reccdf <- function(n, alpha){
  randp <- t(replicate(n, rgamma(length(alpha), alpha)))
  randp <- randp/rowSums(randp)
  randp <- apply(randp,1,p2ccdf)
  return(randp)
}

#' @export confint_dic
confint_dic <- function(out_em, prob=0.95){
  ind=with(out_em, rmlast(event)>0)
  z <- qnorm(1-0.5*(1-prob))
  b <- rmlast(with(out_em, rev(exp(z*sqrt(cumsum(rev(variance[ind])))/log(ccdf[ind])))))
  value=out_em$value[ind]
  ccdf=rmlast(out_em$ccdf)[ind]
  return(data.frame(value=value, ccdf=ccdf,
                    lower=ccdf^(1/b), upper=ccdf^b))
}