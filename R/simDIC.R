#' @export simDIC
simDIC <- function(y, at=NULL, WS=NULL, WE=NULL){
  N <- length(y)
  if(is.null(at)){
    at <- cumsum(rexp(N))
  }
  if(is.null(WS)){
    WS <- rexp(N)
  }
  if(is.null(WE)){
    WE <- rexp(N)
  }
  S <- at + y
  r_e <- runif(N)
  r_s <- runif(N)
  LS <- S - r_s*WS
  RS <- S + (1-r_s)*WS
  LE <- at - r_e*WE
  RE <- at + (1-r_e)*WE
  return(data.frame(LE=pmax(LE,0), RE=RE,
                    LS=pmax(LS,LE), RS=RS,
                    S=S, E=at))
}

