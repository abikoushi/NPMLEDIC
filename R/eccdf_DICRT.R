p2ccdf <- function(p){
  rev(cumsum(rev(p)))
}

#' @export eccdf_DIC
eccdf_DIC <- function(LE, RE, LS, RS, ctype, breaks=NULL, iter=1000L){
  if(is.null(breaks)){
    breaks <- sort(unique(c(RS-RE, RS-LE, LS-RE, LS-LE)))
    breaks <- breaks[breaks>=0]
  }
  n <- length(LE)
  if(length(ctype)==1){
    ctype = rep(ctype, n)
  }
  p <- ep_DIC_em(LE, RE, LS, RS, ctype, breaks, iter)
  ind <- p>0
  return(data.frame(value = breaks[ind], prob = p2ccdf(p[ind])))
}

#' @export eccdf_DICT
eccdf_DICT <- function(LE, RE, LS, RS, ctype,
                       tmax=Inf,
                       breaks=NULL, iter=1000L){
  if(is.null(breaks)){
    breaks <- sort(unique(c(RS-RE, RS-LE, LS-RE, LS-LE)))
    breaks <- breaks[breaks>=0]
  }
  n <- length(EL)
  if(length(ctype)==1){
    ctype = rep(ctype, n)
  }
  if(length(tmax)==1){
    tmax = rep(tmax, n)
  }
  p <- ep_DICT_em(LE, RE, LS, RS, tmax, ctype, breaks, iter)
  ind <- p>0
  return(data.frame(value = breaks[ind], prob = p2ccdf(p[ind])))
}
