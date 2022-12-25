#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// void a_up1(arma::mat & a, const arma::vec & p, const arma::umat & ind) {
//     int n = ind.n_rows;
//     a.fill(0);
//     for(int i=0;i<n;i++){
//     arma::vec q = arma::zeros<arma::vec>(ind(i,1)-ind(i,0)+1);
//         for(int j = arma::as_scalar(ind(i,0));j<arma::as_scalar(ind(i,1));j++){
//             q.row(j-arma::as_scalar(ind(i,0))) += p.row(j);
//         }
//     q /= sum(q);
//     a.rows(ind(i,0), ind(i,1)) += q;
//     }
// }


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

// void a_up2(arma::mat & a, const arma::vec & p, const arma::umat & Lind, const arma::umat & Rind) {
//     int n = Lind.n_rows;
//     a.fill(0);
//     for(int i=0;i<n;i++){
//         arma::vec q = arma::zeros<arma::vec>(Rind(i,1)-Lind(i,0)+1);
//         if(Lind(i,0) < Rind(i,0)){
//             for(int k=arma::as_scalar(Rind(i,0));k<arma::as_scalar(Rind(i,1));k++){
//                 for(int j=arma::as_scalar(Lind(i,0));j<k;j++){
//                     q.row(j-arma::as_scalar(Lind(i,0))) += p.row(j);
//                     }
//             }
//         double den = sum(q);
//         if(den>0){
//             q /= den;
//             a.rows(Lind(i,0), Rind(i,1)) += q;   
//         }
//     }
//     }
// }

arma::vec a_up3(const arma::vec & p, const int & Lind1, const int & Lind2, const int & Rind1, const int & Rind2) {
    arma::vec a = arma::zeros<arma::vec>(p.n_rows);
    arma::vec q = arma::zeros<arma::vec>(Rind2-Lind1+1);
    if(Lind1 < Rind1){
        for(int k=Lind1; k<=Rind1; k++){
            for(int j=Lind1; j<=k; j++){
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
        //Rprintf("%d: ",i);
        // Rprintf("%d,",Lind(i,0));
        // Rprintf("%d\n",Rind(i,0));
        if(ctype[i] == 0){
            a(Lind(i,0)) += 1.0;
        }else if(ctype[i] == 1){
            a += a_up1(p, Lind(i,0), Rind(i,1));
        }else if(ctype[i] == 2){
            a += a_up1(p, Rind(i,0), Rind(i,1));
        }else if(ctype[i] == 3){
            a += a_up3(p, Lind(i,0), Lind(i,1), Rind(i,0), Rind(i,1));
        }
        //Rprintf("\n");
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

void p_up(arma::vec & p, const arma::vec & d) {
    arma::vec num = d;
    p = num/sum(num);
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

// arma::vec ep_ICRT_em(const arma::vec & EL,
//                       const arma::vec & ER,
//                       const arma::vec & S,
//                       const double & tmax,
//                       const arma::vec & breaks,
//                       const int & iter) {
//     int m = breaks.n_rows;
//     arma::vec prob = arma::ones<arma::vec>(m)/m;
//     arma::umat aind = acount(S-ER, S-EL, breaks);
//     arma::umat bind = bcount(tmax-S, breaks);
//     arma::vec d = arma::zeros<arma::vec>(m);
//     arma::vec B = arma::zeros<arma::vec>(m);
//     for(int it=0; it<iter; it++){
//         a_up1(d, prob, aind);
//         b_up(B, d, prob, bind.col(0), bind.col(1));
//         p_up(prob, B, d);
//     }
//     return prob;
// }

// arma::vec ep_DICRT_em(const arma::vec & EL,
//                      const arma::vec & ER,
//                      const arma::vec & SL,
//                      const arma::vec & SR,
//                      const double & tmax,
//                      const arma::vec & breaks,
//                      const int & iter) {
//     int m = breaks.n_rows;
//     arma::vec prob = arma::ones<arma::vec>(m)/m;
//     arma::umat aind_R = acount(SL-EL, SR-EL, breaks);
//     arma::umat aind_L = acount(SL-ER, SR-ER, breaks);
//     arma::umat bind = bcount(tmax-SR, breaks);
//     arma::vec A = arma::zeros<arma::vec>(m);
//     arma::vec B = arma::zeros<arma::vec>(m);
//     for(int it=0; it<iter; it++){
//         a_up2(A, prob, aind_L, aind_R);
//         b_up(B, A, prob, bind.col(0), bind.col(1));
//         p_up(prob, B, A);
//     }
//     return prob;
// }


// [[Rcpp::export]]
arma::vec ep_DIC_em(const arma::vec & EL,
                    const arma::vec & ER,
                    const arma::vec & SL,
                    const arma::vec & SR,
                    const arma::uvec & ctype,
                    const arma::vec & breaks,
                    const int & iter) {
    int m = breaks.n_rows;
    arma::vec prob = arma::ones<arma::vec>(m)/m;
    //     arma::umat aind = acount(S-ER, S-EL, breaks);
    //     arma::umat aind_R = acount(SL-EL, SR-EL, breaks);
    //     arma::umat aind_L = acount(SL-ER, SR-ER, breaks);
    arma::umat aind_R = acount(SL-EL, SR-EL, breaks);
    arma::umat aind_L = acount(SL-ER, SR-ER, breaks);
    arma::vec alpha = arma::zeros<arma::vec>(m);
    for(int it=0; it<iter; it++){
        a_up(alpha, prob, aind_L, aind_R, ctype);
        p_up(prob, alpha);
    }
    return prob;
}

// [[Rcpp::export]]
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
