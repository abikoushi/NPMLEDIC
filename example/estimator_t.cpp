#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

void p_up(arma::vec & p, const arma::vec & b, const arma::vec & d) {
  arma::vec num = d + b;
  p = num/sum(num);
}

void b_up(arma::vec & b, const arma::vec & d, const arma::vec & p, const arma::uvec & j1, const arma::uvec & j2) {
  int m = b.n_elem;
  int n = j1.n_elem;
  b.fill(0);
  for(int i=0;i<n;i++){
    double q = sum(p.rows(j1[i], j2[i]));
    if(q > 0.0){
      double bj = (1-q)/q;
      for(int j=0; j<j1[i]; j++){
        b.row(j) += bj*d.row(j);   
      }
      for(int j=j2[i]; j<m; j++){
        b.row(j) += bj*d.row(j);   
      }
    }
  }
}

arma::umat bcount(const arma::vec & L, const arma::vec & U, const arma::vec & breaks){
  int n = U.n_rows;
  int m = breaks.n_rows;
  arma::umat indexmat = arma::zeros<arma::umat>(n, 2);
  for(int i=0; i<n; i++){
    bool start = true;
    for(int j=0; j<m; j++){
      bool beta = arma::as_scalar(L.row(i) <= breaks.row(j) && breaks.row(j) <= U.row(i));
      if(start){
        if(beta){
          indexmat(i,0) = j;
          indexmat(i,1) = j;
          start = false;
        }
      }else{
        if(beta){
          indexmat(i,1) = j;
        }
      }
    }
  }
  return indexmat;
}


arma::umat bcount(const arma::vec & U, const arma::vec & breaks){
  int n = U.n_rows;
  int m = breaks.n_rows;
  arma::umat indexmat = arma::zeros<arma::umat>(n, 2);
  for(int i=0; i<n; i++){
    bool start = true;
    for(int j=0; j<m; j++){
      bool beta = arma::as_scalar(breaks.row(j) <= U.row(i));
      if(start){
        if(beta){
          indexmat(i,0) = j;
          indexmat(i,1) = j;
          start = false;
        }
      }else{
        if(beta){
          indexmat(i,1) = j;
        }
      }
    }
  }
  return indexmat;
}


arma::vec ep_DICT_em(const arma::vec & EL,
                    const arma::vec & ER,
                    const arma::vec & SL,
                    const arma::vec & SR,
                    const arma::vec & tmax,
                    const arma::uvec & ctype,
                    const arma::vec & breaks,
                    const int & iter) {
    int m = breaks.n_rows;
    arma::vec prob = arma::ones<arma::vec>(m)/m;
    //  arma::umat aind = acount(S-ER, S-EL, breaks);
    //  arma::umat aind_R = acount(SL-EL, SR-EL, breaks);
    //  arma::umat aind_L = acount(SL-ER, SR-ER, breaks);
    arma::umat aind_R = acount(SL-EL, SR-EL, breaks);
    arma::umat aind_L = acount(SL-ER, SR-ER, breaks);
    arma::umat bind_R = bcount(tmax-EL, breaks);
    arma::umat bind_L = bcount(tmax-ER, breaks);
    arma::vec Alpha = arma::zeros<arma::vec>(m);
    arma::vec Beta = arma::zeros<arma::vec>(m);
    for(int it=0; it<iter; it++){
        a_up(Alpha, prob, aind_L, aind_R, ctype);
        b_up(Beta, Alpha, prob, bind_L, bind_R);
        p_up(prob, Alpha);
    }
    return prob;
}
