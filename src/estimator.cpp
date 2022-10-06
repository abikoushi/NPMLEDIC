// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


void d_up(arma::mat & d, const arma::vec & p, const arma::uvec & j1, const arma::uvec & j2) {
    //int m = d.n_elem;
    int n = j1.n_elem;
    d.fill(0);
    for(int i=0;i<n;i++){
        Rprintf("%d ",i);
        if(j1[i]>=0 && j2[i]>=0){
            arma::vec q = p.rows(j1[i], j2[i]);
            q /= sum(q);
            d.rows(j1[i],j2[i]) += q;
        }
    }
}

void b_up(arma::vec & b, const arma::vec & d, const arma::vec & p, const arma::uvec & j1, const arma::uvec & j2) {
    int m = b.n_elem;
    int n = j1.n_elem;
    b.fill(0);
    for(int i=0;i<n;i++){
        if(j1[i]>=0 && j2[i]>=0){
            double q = sum(p.rows(j1[i],j2[i]));
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
}


void p_up(arma::vec & p, const arma::vec & b, const arma::vec & d) {
    arma::vec num = d + b;
    p = num/sum(num);
}

arma::umat acount(const arma::vec & L, const arma::vec & R, const arma::vec & breaks){
    int n = L.n_rows;
    int m = breaks.n_rows;
    arma::umat indexmat = arma::zeros<arma::umat>(n,2);
    for(int i=0; i<n; i++){
    bool st = true;
    for(int j=0; j<m-1; j++){
        bool alpha = arma::as_scalar(L.row(i) <= breaks.row(j) && breaks.row(j+1) <= R.row(i));
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

        
arma::umat bcount(const arma::vec & U, const arma::vec & breaks){
    int n = U.n_rows;
    int m = breaks.n_rows;
    arma::umat indexmat = arma::zeros<arma::umat>(n, 2);
    for(int i=0; i<n; i++){
        bool st = true;
        for(int j=0; j<m; j++){
            bool beta = arma::as_scalar(breaks.row(j) <= U.row(i));
            if(st){
                if(beta){
                    indexmat(i,0) = j;
                    indexmat(i,1) = j;
                    st = false;
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

// arma::vec getS(const arma::vec & SL, const arma::vec & SR, const double & midp){
//     return SL+midp*(SR-SL);
// }


// [[Rcpp::export]]
arma::vec ep_ICRT_em(const arma::vec & EL,
                      const arma::vec & ER,
                      const arma::vec & S,
                      const double & tmax,
                      const arma::vec & breaks,
                      const int & iter) {
    int m = breaks.n_rows;
    arma::vec prob = arma::ones<arma::vec>(m)/m;
    arma::umat aind = acount(S-ER, S-EL, breaks);
    arma::umat bind = bcount(tmax-S, breaks);
    arma::vec d = arma::ones<arma::vec>(m);
    arma::vec B = arma::ones<arma::vec>(m);
    aind.print();
    for(int it=0; it<iter; it++){
        d_up(d, prob, aind.col(0), aind.col(1));
        b_up(B, d, prob, bind.col(0), bind.col(1));
        p_up(prob, B, d);
    }
    return prob;
}
