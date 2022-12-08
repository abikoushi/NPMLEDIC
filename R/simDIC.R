#' @export simDIC
simDIC <- function(y){
  N <- length(y)
  at <- cumsum(rexp(N))
  Sw <- rexp(N)
  Ew <- rexp(N)
  S <- at + y
  r_e <- runif(N)
  r_s <- runif(N)
  LS <- S - r_s*Sw
  RS <- S + (1-r_s)*Sw
  LE <- at - r_e*Ew
  RE <- at + (1-r_e)*Ew
  return(data.frame(LE=pmax(LE,0), RE=RE,
                    LS=pmax(LS,LE), RS=RS,
                    S=S, E=at))
}

