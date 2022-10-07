#' @export eccdf_DICRT
eccdf_DICRT <- function(EL, ER, SL, SR, tmax=Inf, breaks=NULL, iter=1000L){
  if(is.null(breaks)){
    breaks <- sort(unique(c(SR-ER, SR-EL, SL-ER, SL-EL)))
    breaks <- breaks[breaks>=0]
  }
  p <- ep_DICRT_em(EL, ER, SL, SR, tmax, breaks, iter) 
  return(data.frame(value = breaks, prob = rev(cumsum(rev(p)))))
}

#' @export eccdf_ICRT
eccdf_ICRT <- function(EL, ER, S, tmax=Inf, breaks=NULL, iter=1000L){
  if(is.null(breaks)){
    breaks <- sort(unique(c(S-ER, S-EL)))
    breaks <- breaks[breaks>=0]
  }
  p <- ep_ICRT_em(EL, ER, S, tmax, breaks, iter) 
  return(data.frame(value = breaks, prob = rev(cumsum(rev(p)))))
}
