#' @export eccdf_DICRT
eccdf_DICRT <- function(EL, ER, SL, SR, tmax=Inf, breaks=NULL, midp=0.5, iter=1000L){
  S <- SL + midp*(SR-SL)
  if(is.null(breaks)){
    breaks <- sort(unique(c(S-ER, S-EL)))
  }
  p <- ep_ICRT_em(EL, ER, S, tmax, breaks, iter) 
  return(data.frame(value = breaks, ccdf = 1-cumsum(p)))
}
