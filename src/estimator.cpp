#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

void a_up1(arma::mat & a, const arma::vec & p, const arma::umat & ind) {
    int n = ind.n_rows;
    a.fill(0);
    for(int i=0;i<n;i++){
    arma::vec q = arma::zeros<arma::vec>(ind(i,1)-ind(i,0)+1);
        for(int j = arma::as_scalar(ind(i,0));j<arma::as_scalar(ind(i,1));j++){
            q.row(j-arma::as_scalar(ind(i,0))) += p.row(j);
        }
    q /= sum(q);
    a.rows(ind(i,0), ind(i,1)) += q;
    }
}

void a_up2(arma::mat & a, const arma::vec & p, const arma::umat & Lind, const arma::umat & Rind) {
    int n = Lind.n_rows;
    a.fill(0);
    for(int i=0;i<n;i++){
        arma::vec q = arma::zeros<arma::vec>(Rind(i,1)-Lind(i,0)+1);
        for(int k=arma::as_scalar(Rind(i,0));k<arma::as_scalar(Rind(i,1));k++){
        for(int j=arma::as_scalar(Lind(i,0));j<k;j++){
            q.row(j-arma::as_scalar(Lind(i,0))) += p.row(j);
        }
        }
        double den = sum(q);
        if(den>0){
            q /= den;
            a.rows(Lind(i,0), Rind(i,1)) += q;   
        }
    }
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
    arma::vec d = arma::zeros<arma::vec>(m);
    arma::vec B = arma::zeros<arma::vec>(m);
    for(int it=0; it<iter; it++){
        a_up1(d, prob, aind);
        b_up(B, d, prob, bind.col(0), bind.col(1));
        p_up(prob, B, d);
    }
    return prob;
}

// [[Rcpp::export]]
arma::vec ep_DICRT_em(const arma::vec & EL,
                     const arma::vec & ER,
                     const arma::vec & SL,
                     const arma::vec & SR,
                     const double & tmax,
                     const arma::vec & breaks,
                     const int & iter) {
    int m = breaks.n_rows;
    arma::vec prob = arma::ones<arma::vec>(m)/m;
    arma::umat aind_R = acount(SL-EL, SR-EL, breaks);
    arma::umat aind_L = acount(SL-ER, SR-ER, breaks);
    arma::umat bind = bcount(tmax-SR, breaks);
    arma::vec A = arma::zeros<arma::vec>(m);
    arma::vec B = arma::zeros<arma::vec>(m);
    for(int it=0; it<iter; it++){
        a_up2(A, prob, aind_L, aind_R);
        b_up(B, A, prob, bind.col(0), bind.col(1));
        p_up(prob, B, A);
    }
    return prob;
}
