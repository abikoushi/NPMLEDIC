#' @export simDIC
simDIC <- function(inc, Ew, Sw){
  N <- length(inc)
  r_e <- runif(N)
  r_s <- runif(N)
  S_L <- inc - r_e*Sw
  S_R <- inc + (1-r_e)*Sw
  E_R <- (1-r_s)*Ew
  data.frame(E_L=0, E_R=pmin(E_R,S_R),
             S_L=pmax(S_L,0), S_R=S_R)
}