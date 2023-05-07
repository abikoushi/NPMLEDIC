
p2ccdf <- function(p){
  rev(cumsum(rev(p)))
}

rmlast <- function(x)x[-length(x)]

setbreaks <- function(LE, RE, LS, RS){
  #LS-RE #shortest case
  #RS-LE #longest case
  breaks <- sort(unique(c(0, LS-RE, RS-LE)))
  return(breaks[breaks>=0])
}


#' @export eccdf_dic_em
eccdf_dic_em <- function(LE, RE, LS, RS, ctype,
                         alpha0 = 1, maxit=1000L, tol=1e-5){
  breaks <- setbreaks(LE, RE, LS, RS)
  n <- length(LE)
  if(length(ctype)==1L){
    ctype = rep(ctype, n)
  }
  res <- ep_DIC_em(LE, RE, LS, RS, ctype, breaks, alpha0, maxit, tol)
  res$value <- breaks
  res$ccdf  <- with(res, p2ccdf(prob))
  m <- length(res$prob)
  I <- matrix(1,m,m)
  diag(I) <- 2/res$prob-1
  res$I <- sum(res$event)*I
  return(res)
}

#' @export confint_dic
confint_dic <- function(out_em, prob=0.95, scale="linear"){
  ind=with(out, event>0)
  variance <- diag(solve(out$I[ind,ind]))
  z <- qnorm(1-0.5*(1-probs))
  value = out_em$value[ind]
  ccdf = out$ccdf[ind]
  if(scale=="linear"){
    b <- sapply(z, function(z0)rev(z0*sqrt(cumsum(rev(variance)))))
    lower = as.matrix(ifelse(ccdf-b<0, 0, ccdf-b))
    upper = as.matrix(ifelse(ccdf+b>1, 1, ccdf+b))
  }
  if(scale=="loglog"){
    b <- sapply(z, function(z0)rev(z0*sqrt(cumsum(rev(variance)))/log(ccdf)))
    lower = as.matrix(ccdf^exp(-b))
    upper = as.matrix(ccdf^exp(b))
  }
  res <- lapply(1:length(probs),
                function(i)data.frame(value=value,
                                      lower = lower[,i],
                                      upper = upper[,i],
                                      level = probs[i]))
  res <- do.call("rbind",res)
  return(list(point=data.frame(value=value,ccdf=ccdf),
              confint=res))
}

#' @export eccdf_dic_gibbs
eccdf_dic_gibbs <- function(LE, RE, LS, RS, ctype, alpha0 = 0, iter=2000L){
  breaks <- setbreaks(LE, RE, LS, RS)
  n <- length(LE)
  if(length(ctype)==1){
    ctype = rep(ctype, n)
  }
  res <- ep_DIC_gibbs(LE, RE, LS, RS, ctype, breaks, alpha0, iter)    
  res$value <- breaks
  res$ccdf  <- with(res, apply(prob, 1, p2ccdf))
  return(res)
}


#' @export eccdf_dic_vb
eccdf_dic_vb <- function(LE, RE, LS, RS, ctype,
                         alpha0 = 1, iter = 1000L){
  breaks <- setbreaks(LE, RE, LS, RS)
  
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
  randp <- apply(randp, 1, p2ccdf)
  return(randp)
}

