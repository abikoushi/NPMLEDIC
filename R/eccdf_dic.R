
p2ccdf <- function(p){
  rev(cumsum(rev(p)))
}

rmlast <- function(x)x[-length(x)]

setbreaks <- function(LE, RE, LS, RS){
  #LS-RE #shortest case
  #RS-LE #longest case
  # breaks <- sort(unique(c(0, LS-RE, RS-LE))) #old
  breaks <- sort(unique(c(LS-RE, RS-LE)))
  return(breaks[breaks>=0])
}


acount <- function(LS,RS,E,breaks){
  #the index is 0-origin in cpp (-1)   
  cbind(match(LS-E,breaks), match(RS-E,breaks))-1L
}

#' @export jointecdf_dic_em
jointecdf_dic_em <- function(LE, RE, LS, RS, maxit=1000L, tol=1e-5){
  breaks <- sort(unique(c(LE, RE, LS, RS)))
  LE <- match(LE,breaks)
  RE <- match(RE,breaks)
  LS <- match(LS,breaks)
  RS <- match(RS,breaks)
  res <- joint_DIC_em(LE-1L, RE-1L, LS-1L, RS-1L, length(breaks), maxit, tol)
  res$breaks <- breaks
  return(res)
}

#' @export eccdf_dic_em
eccdf_dic_em <- function(LE, RE, LS, RS,
                         alpha0=0,
                         maxit=1000L, tol=1e-5){
  breaks <- setbreaks(LE, RE, LS, RS)
  alpha0 <- diff(c(0,breaks))*alpha0 #prop. to interval-length
  n <- length(LE)
  if(length(ctype)==1L){
    ctype = rep(ctype, n)
  }
  aind_L <- acount(LS,RS,RE,breaks)
  aind_R <- acount(LS,RS,LE,breaks)
  # aind_L <- cbind(match(LS-RE,breaks), match(RS-RE,breaks))-1L
  # aind_R <- cbind(match(LS-LE,breaks), match(RS-LE,breaks))-1L
  res <- ep_DIC_em(aind_L, aind_R, ctype, alpha0, maxit, tol)
  #res <- ep_DIC_em(LE, RE, LS, RS, ctype, breaks, alpha0, maxit, tol)
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
  breaks <- setbreaks(LE, RE, LS, RS)
  alpha0 <- diff(c(0,breaks))*beta #prop. to interval-length
  n <- length(LE)
  if(length(ctype)==1L){
    ctype = rep(ctype, n)
  }
  aind_L <- acount(LS,RS,RE,breaks)
  aind_R <- acount(LS,RS,LE,breaks)
  n <- length(LE)
  if(length(ctype)==1L){
    ctype = rep(ctype, n)
  }
  res <- ep_DIC_vb(LE, RE, LS, RS, ctype, alpha0, maxit, tol)
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
eccdf_dic_gibbs <- function(LE, RE, LS, RS, ctype, beta = 1, iter=2000L){
  breaks <- setbreaks(LE, RE, LS, RS)
  breaks <- setbreaks(LE, RE, LS, RS)
  alpha0 <- diff(c(0,breaks))*beta #prop. to interval-length
  n <- length(LE)
  if(length(ctype)==1L){
    ctype = rep(ctype, n)
  }
  aind_L <- acount(LS,RS,RE,breaks)
  aind_R <- acount(LS,RS,LE,breaks)
  n <- length(LE)
  res <- ep_DIC_gibbs(LE, RE, LS, RS, ctype, breaks, alpha0, iter)    
  res$value <- breaks
  res$ccdf  <- with(res, apply(prob, 1, p2ccdf))
  return(res)
}
