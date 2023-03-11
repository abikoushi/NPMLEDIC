#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double sumxlogy(const arma::vec & x, const arma::vec & y){
    double out=0;
    for(int i=0; i<x.n_rows; i++){
        if(x[i]>0){
            out += x[i]*log(y[i]);
        }
    }
    return out;
}

arma::vec vec_digamma(arma::vec a){
    int K = a.n_rows;
    arma::vec out(K);
    for(int k=0;k<K;k++){
        out[k] = R::digamma(a[k]);
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


double p_up(arma::vec & p, const arma::vec & d) {
    arma::vec num = d;
    p = num/sum(num);
    return sumxlogy(d, p);
}

double elp_up(arma::vec & elp, arma::vec & alpha, const arma::vec & d, const double & alpha0) {
    alpha = d + alpha0;
    arma::vec logp = vec_digamma(alpha)-R::digamma(sum(alpha));
    elp = exp(logp);
    return sum(d%logp);
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

// [[Rcpp::export]]
List ep_DIC_em(const arma::vec & EL,
                    const arma::vec & ER,
                    const arma::vec & SL,
                    const arma::vec & SR,
                    const arma::uvec & ctype,
                    const arma::vec & breaks,
                    const int & iter) {
    int m = breaks.n_rows;
    arma::vec prob = arma::ones<arma::vec>(m+1)/(m+1);
    arma::umat aind_R = acount(SL-EL, SR-EL, breaks);
    arma::umat aind_L = acount(SL-ER, SR-ER, breaks);
    arma::vec d = arma::zeros<arma::vec>(m+1);
    arma::vec lp = arma::zeros<arma::vec>(iter);
    for(int it=0; it<iter; it++){
        a_up(d, prob, aind_L, aind_R, ctype);
        lp.row(it) = p_up(prob, d);
    }
    return List::create(_["prob"]=prob, _["event"]=d, _["lp"]=lp);
}

// [[Rcpp::export]]
List ep_DIC_vb(const arma::vec & EL,
                    const arma::vec & ER,
                    const arma::vec & SL,
                    const arma::vec & SR,
                    const arma::uvec & ctype,
                    const arma::vec & breaks,
                    const double & alpha0,
                    const int & iter) {
    int m = breaks.n_rows;
    arma::vec alpha = arma::ones<arma::vec>(m+1);
    arma::vec prob = alpha/(m+1);
    arma::umat aind_R = acount(SL-EL, SR-EL, breaks);
    arma::umat aind_L = acount(SL-ER, SR-ER, breaks);
    arma::vec d = arma::zeros<arma::vec>(m+1);
    arma::vec lp = arma::zeros<arma::vec>(iter);
    for(int it=0; it<iter; it++){
        a_up(d, prob, aind_L, aind_R, ctype);
        lp.row(it) = elp_up(prob, alpha, d, alpha0);
    }
    return List::create(_["alpha"]=alpha, _["event"]=d, _["lp"]=lp);
}
