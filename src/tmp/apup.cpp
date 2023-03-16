#include <RcppArmadillo.h>
#include "apup.h"

double sumxlogy(const arma::vec & x, const arma::vec & y){
  double out=0;
  for(int i=0; i<x.n_rows; i++){
    if(x[i]>0){
      out += x[i]*log(y[i]);
    }
  }
  return out;
}

//interval censored E
arma::vec a_up1(const arma::vec & p, const int & ind1, const int & ind2) {
  arma::vec a = arma::zeros<arma::vec>(p.n_rows);
  if(ind1 < ind2){
    arma::vec q = arma::zeros<arma::vec>(ind2-ind1+1);
    for(int j = ind1;j<=ind2;j++){
      //Rprintf("%d ", j);
      q.row(j-ind1) += p.row(j);
    }
    q /= sum(q);
    a.rows(ind1, ind2) += q;
  }
  return a;
}

//interval censored both
arma::vec a_up3(const arma::vec & p, const int & Lind1, const int & Lind2, const int & Rind1, const int & Rind2) {
  arma::vec a = arma::zeros<arma::vec>(p.n_rows);
  arma::vec q = arma::zeros<arma::vec>(Rind2-Lind1+1);
  if(Lind1 < Rind1){
    for(int k=Lind1; k<=Rind1; k++){// <= ?
      for(int j=Lind1; j<=k; j++){
        //Rprintf("%d ", j);
        q.row(j-Lind1) += p.row(j);
      }
    }
    double den = sum(q);
    if(den>0){
      q /= den;
      a.rows(Lind1, Rind2) += q;
    }
  }
  return a;
}

arma::umat acount(const arma::vec & L, const arma::vec & R, const arma::vec & breaks){
  int n = L.n_rows;
  int m = breaks.n_rows;
  arma::umat indexmat = arma::zeros<arma::umat>(n,2);
  for(int i=0; i<n; i++){
    bool st = true;
    for(int j=0; j<m; j++){
      bool alpha = arma::as_scalar(L.row(i) <= breaks.row(j) && breaks.row(j) <= R.row(i));
      if(st){
        if(alpha){
          indexmat(i,0) = j;
          indexmat(i,1) = j;
          st = false;
        }
      }else{
        if(alpha){
          indexmat(i,1) = j;
        }
      }
    }
  }
  return indexmat;
}

double p_up(arma::vec & p, const arma::vec & d, const double & alpha0) {
  arma::vec num = d + alpha0;
  p = num/sum(num);
  return sumxlogy(d, p);
}
