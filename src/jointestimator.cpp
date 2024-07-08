#include <RcppArmadillo.h>
#include "util.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec expmeanlog(arma::vec x){
  double sumx = sum(x);
  arma::vec out = arma::zeros<arma::vec>(x.n_rows);
  //out.fill(0);
  for(int i=0; i<x.n_rows; i++){
    if(x(i)>0){
      out.row(i) = exp(R::digamma(x(i)) - R::digamma(sumx)); 
    }
  }
  return out;
}

arma::mat a_upsub(const arma::vec & h, const arma::vec & lam,
                const int & lerank, const int & rerank,
                const int & lsrank, const int & rsrank) {
  arma::mat a = arma::zeros<arma::mat>(h.n_rows, h.n_rows);
  for(int i=lerank; i<rsrank; i++){
    //Rprintf("i %d\n", i);
    int sta = std::max(lsrank,i);
    for(int j=sta; j<rsrank; j++){
      //Rprintf("j %d\n", j);
        a.col(j).row(i) += h.row(j-i)*lam.row(i);
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

double par_up(arma::vec & alpha_h, arma::vec & beta_l, arma::vec & h, arma::vec & lam, const arma::mat & A) {
  int m = h.n_rows;
  alpha_h.fill(0);
  beta_l.fill(0);
  for(int e=0; e<m; e++){
    for(int s=e; s<m; s++){
      alpha_h.row(s-e) += A(e,s);
      beta_l.row(e) += A(e,s);
    }
  }
  h = expmeanlog(alpha_h);
  lam = expmeanlog(beta_l);
  return sumxlogy(alpha_h, h)+sumxlogy(beta_l, lam);
}

// [[Rcpp::export]]
List joint_DIC_em(const arma::uvec & EL,
               const arma::uvec & ER,
               const arma::uvec & SL,
               const arma::uvec & SR,
               const int & m,
               const int & maxit, const double & tol) {
  arma::vec alpha_h = arma::ones<arma::vec>(m);
  arma::vec beta_l = arma::ones<arma::vec>(m);
  arma::mat A = arma::zeros<arma::mat>(m,m);
  arma::vec lp = arma::zeros<arma::vec>(maxit);
  arma::vec h = arma::ones<arma::vec>(m)/m;
  arma::vec lam = arma::ones<arma::vec>(m)/m;
  for(int it=1; it<maxit; it++){
    Aup(A, h, lam, EL, ER, SL, SR);
    lp.row(it) = par_up(alpha_h, beta_l, h, lam, A);
    if(abs(arma::as_scalar(lp.row(it)-lp.row(it-1))) < tol){
      lp = lp.rows(1,it);
      break;
    }
  }
  return List::create(_["alpha"]=alpha_h,
                      _["beta"]=beta_l,
                      _["event"]=A,
                      _["lp"]=lp);
}
