p2ccdf <- function(p){
  rev(cumsum(rev(p)))
}

#' @export eccdf_DIC
eccdf_DIC <- function(EL, ER, SL, SR, ctype, breaks=NULL, iter=1000L){
  if(is.null(breaks)){
    breaks <- sort(unique(c(SR-ER, SR-EL, SL-ER, SL-EL)))
    breaks <- breaks[breaks>=0]
  }
  p <- ep_DIC_em(EL, ER, SL, SR, ctype, breaks, iter)
  ind <- p>0
  return(data.frame(value = breaks[ind], prob = p2ccdf(p[ind])))
}


#' eccdf_DICRT <- function(EL, ER, SL, SR, tmax=Inf, breaks=NULL, iter=1000L){
#'   if(is.null(breaks)){
#'     breaks <- sort(unique(c(SR-ER, SR-EL, SL-ER, SL-EL)))
#'     breaks <- breaks[breaks>=0]
#'   }
#'   p <- ep_DICRT_em(EL, ER, SL, SR, tmax, breaks, iter) 
#'   return(data.frame(value = breaks, prob = p2ccdf(p)))
#' }
