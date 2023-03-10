
p2ccdf <- function(p){
  rev(cumsum(rev(p)))
}

#' @export eccdf_dic_em
eccdf_dic_em <- function(LE, RE, LS, RS, ctype, iter=1000L){
  breaks <- sort(unique(c(RS-RE, RS-LE, LS-RE, LS-LE)))
  breaks <- breaks[breaks>=0]
  
  n <- length(LE)
  if(length(ctype)==1){
    ctype = rep(ctype, n)
  }
  p <- ep_DIC_em(LE, RE, LS, RS, ctype, breaks, iter)    
  ind <- p>0
  return(data.frame(value = breaks[ind], prob = p2ccdf(p[ind])))
}

#' @export eccdf_dic_vb
eccdf_dic_vb <- function(LE, RE, LS, RS, ctype, alpha0=0.5, iter=1000L){
  breaks <- sort(unique(c(RS-RE, RS-LE, LS-RE, LS-LE)))
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
