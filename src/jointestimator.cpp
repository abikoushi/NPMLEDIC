#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double sumxlogy2(const arma::vec & x, const arma::vec & y){
  double out=0;
  for(int i=0; i<x.n_rows; i++){
    if(x[i]>0){
      out += x[i]*log(y[i]);
    }
  }
  return out;
}

arma::mat a_upsub(const arma::vec & h, const arma::vec & lam,
                const int & lerank, const int & rerank,
                const int & lsrank, const int & rsrank) {
  arma::mat a = arma::zeros<arma::mat>(h.n_rows, h.n_rows);
  for(int i=lerank; i<=rsrank; i++){
    for(int j=lsrank; j<=rsrank; j++){
      if(i < j){
        a.col(j).row(i) += h.row(j-i)*lam.row(i);
      }
    }
  }
  a /= arma::accu(a);
  return a;
}

void Aup(arma::mat & A, const arma::vec & h, const arma::vec & lam,
          const arma::uvec & lerank, const arma::uvec & rerank,
          const arma::uvec & lsrank, const arma::uvec & rsrank) {
  int n = lsrank.n_rows;
  A.fill(0.0);
  for(int i=0;i<n;i++){
    A += a_upsub(h, lam, lerank[i], rerank[i], lsrank[i], rsrank[i]);
  }
}

double par_up(arma::vec & h, arma::vec & lam, const arma::mat & A) {
  int m = h.n_rows;
  arma::vec num_h = arma::zeros<arma::vec>(m);
  arma::vec num_l = arma::zeros<arma::vec>(m);
  for(int e=0; e<m; e++){
    for(int s=e+1; s<m; s++){
      num_h.row(s-e) += A(e,s);
      num_l.row(e) += A(e,s);
    }
  }
  h = num_h/sum(num_h);
  lam = num_l/sum(num_l);
  return sumxlogy2(num_h, h)+sumxlogy2(num_l, lam);
}

// [[Rcpp::export]]
List joint_DIC_em(const arma::uvec & EL,
               const arma::uvec & ER,
               const arma::uvec & SL,
               const arma::uvec & SR,
               const int & m,
               const int & maxit, const double & tol) {
  arma::vec h = arma::ones<arma::vec>(m)/(m);
  arma::vec lam = arma::ones<arma::vec>(m)/(m);
  arma::mat A = arma::zeros<arma::mat>(m,m);
  arma::vec lp = arma::zeros<arma::vec>(maxit);
  for(int it=1; it<maxit; it++){
    Aup(A, h, lam, EL, ER, SL, SR);
    lp.row(it) = par_up(h, lam, A);
    if(abs(arma::as_scalar(lp.row(it)-lp.row(it-1))) < tol){
      lp = lp.rows(1,it);
      break;
    }
  }
  return List::create(_["h"]=h, _["lambda"]=lam, _["event"]=A, _["lp"]=lp);
}
