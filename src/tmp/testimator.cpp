#include "RcppArmadillo.h"
#include "apup.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

void a_up(arma::vec & a, const arma::vec & p, const arma::umat & Lind, const arma::umat & Rind, const arma::uvec & ctype) {
  int n = Lind.n_rows;
  a.fill(0);
  for(int i=0;i<n;i++){
    if(ctype[i] == 0){
      a(Lind(i,0)) += 1.0;
    }else if(ctype[i] == 1){
      a += a_up1(p, Lind(i,0), Rind(i,1));
    }else if(ctype[i] == 2){
      a += a_up1(p, Rind(i,0), Rind(i,1));
    }else if(ctype[i] == 3){
      a += a_up3(p, Lind(i,0), Lind(i,1), Rind(i,0), Rind(i,1));
    }
  }
  a.row(a.n_rows-1) += n-sum(a);
}

///
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

void b_up(arma::vec & b, const arma::vec & d, const arma::vec & p, const arma::uvec & jR) {
  int m = b.n_elem;
  int n = jR.n_elem;
  b.fill(0);
  for(int i=0;i<n;i++){
    double q = sum(p.rows(0, jR[i]));
    if(0 < q & q < 1){
      double bj = (1-q)/q;
      double den = sum(p.rows(jR[i]+1, m-1));
      for(int j=jR[i]+1; j<m; j++){
        b.row(j) += bj*p.row(j)/den;
      }
    }
  }
}


// both side (upper & lower) truncation
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

// right (upper) truncation
arma::umat bcount(const arma::vec & U, const arma::vec & breaks){
  int n = U.n_rows;
  int m = breaks.n_rows;
  arma::umat indexmat = arma::zeros<arma::uvec>(n);
  for(int i=0; i<n; i++){
    bool start = true;
    for(int j=0; j<m; j++){
      bool beta = arma::as_scalar(breaks.row(j) <= U.row(i));
      if(beta){
        indexmat(i) = j;
        break;
      }
    }
  }
  return indexmat;
}

List ep_DICT_em(const arma::vec & EL,
                const arma::vec & ER,
                const arma::vec & SL,
                const arma::vec & SR,
                const arma::vec & tmax,
                const arma::uvec & ctype,
                const arma::vec & breaks,
                const double & alpha0,
                const int & iter) {
  int m = breaks.n_rows;
  arma::vec prob = arma::ones<arma::vec>(m+1)/(m+1);
  arma::umat aind_R = acount(SL-EL, SR-EL, breaks);
  arma::umat aind_L = acount(SL-ER, SR-ER, breaks);
  arma::umat bind_R = bcount(tmax-EL, breaks);
  //arma::umat bind_L = bcount(tmax-ER, breaks);
  arma::vec d = arma::zeros<arma::vec>(m+1);
  arma::vec b = arma::zeros<arma::vec>(m+1);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  for(int it=0; it<iter; it++){
    a_up(d, prob, aind_L, aind_R, ctype);
    b_up(b, d, prob, bind_R);
    lp.row(it) = p_up(prob, b+d, alpha0);
  }
  return List::create(_["prob"]=prob, _["event"]=d, _["event2"]=b, _["lp"]=lp);
}
