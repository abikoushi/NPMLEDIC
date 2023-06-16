
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
                         alpha0 = 1e-8, maxit=1000L, tol=1e-5){
  breaks <- setbreaks(LE, RE, LS, RS)
  n <- length(LE)
  if(length(ctype)==1L){
    ctype = rep(ctype, n)
  }
  res <- ep_DIC_em(LE, RE, LS, RS, ctype, breaks, alpha0, maxit, tol)
  est <- data.frame(value = breaks,
                    event = res$event,
                    prob = res$prob,
                    ccdf = p2ccdf(res$prob),
                    cdf = cumsum(res$prob))
  m <- length(res$prob)
  I <- matrix(1,m,m)
  diag(I) <- 2/res$prob-1
  I <- sum(res$event)*I
  diag(I) <- diag(I) + alpha0/(res$prob^2)
  return(list(estimates=est, I=I, lp=res$lp))
}

#' @export confint_dic
confint_dic <- function(out_em, prob=0.95){
  #scale="linear"
  ind=with(out_em, estimates$event>0)
  variance <- diag(solve(out_em$I[ind,ind]))
  z <- qnorm(1-0.5*(1-prob))
  value = out_em$estimates$value[ind]
  ccdf = out_em$estimates$ccdf[ind]
  #if(scale=="linear"){
    # b <- sapply(z, function(z0)rev(z0*sqrt(cumsum(rev(variance)))))
    b <- sapply(z, function(z0)rev(z0*sqrt(cumsum(variance))))
    lower = as.matrix(ifelse(ccdf-b<0, 0, ccdf-b))
    upper = as.matrix(ifelse(ccdf+b>1, 1, ccdf+b))
  #}
  # if(scale=="loglog"){
  #   b <- sapply(z, function(z0)rev(z0*sqrt(cumsum(rev(variance)))/log(ccdf)))
  #   lower = as.matrix(ccdf^exp(-b))
  #   upper = as.matrix(ccdf^exp(b))
  # }
  res <- lapply(1:length(prob),
                function(i)data.frame(value=value,
                                      ccdf=ccdf,
                                      se=b[,i],
                                      lower = lower[,i],
                                      upper = upper[,i],
                                      level = prob[i]))
  res <- do.call("rbind",res)
  #if(!return_var){
    return(res) 
  # }else{
  #   return(list(CI=res, variance=data.frame(value = value, ccdf=ccdf, se=b)))
  # }
}

######
#VB tools

#' @export eccdf_dic_vb
eccdf_dic_vb <- function(LE, RE, LS, RS, ctype,
                         beta = 1, maxit = 1000L, tol=1e-5){
  breaks <- setbreaks(LE, RE, LS, RS)
  n <- length(LE)
  if(length(ctype)==1L){
    ctype = rep(ctype, n)
  }
  res <- ep_DIC_vb(LE, RE, LS, RS, ctype, breaks, beta, maxit, tol)
  prob <- with(res, alpha/sum(alpha))
  # res$ccdf  <- p2ccdf(prob)
  # res$cdf <- cumsum(prob)
  # res$value <- breaks
  est <- data.frame(value = breaks,
                    alpha = res$alpha,
                    event = res$event,
                    prob = prob,
                    ccdf = p2ccdf(prob),
                    cdf = cumsum(prob))
  
  return(list(estimates=est, lp=res$lp))
}

#' @export reccdf
reccdf <- function(n, alpha){
  randp <- t(replicate(n, rgamma(length(alpha), alpha)))
  randp <- randp/rowSums(randp)
  randp <- apply(randp, 1, p2ccdf)
  return(randp)
}

####
#supplemental

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
