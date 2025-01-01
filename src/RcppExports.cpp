// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ep_DIC_em
List ep_DIC_em(const arma::umat aind_L, const arma::umat aind_R, const arma::uvec& ctype, const arma::vec& alpha0, const int& maxit, const double& tol);
RcppExport SEXP _NPMLEDIC_ep_DIC_em(SEXP aind_LSEXP, SEXP aind_RSEXP, SEXP ctypeSEXP, SEXP alpha0SEXP, SEXP maxitSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat >::type aind_L(aind_LSEXP);
    Rcpp::traits::input_parameter< const arma::umat >::type aind_R(aind_RSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ctype(ctypeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(ep_DIC_em(aind_L, aind_R, ctype, alpha0, maxit, tol));
    return rcpp_result_gen;
END_RCPP
}
// ep_DIC_vb
List ep_DIC_vb(const arma::umat aind_L, const arma::umat aind_R, const arma::uvec& ctype, const arma::vec& alpha0, const int& maxit, const double& tol);
RcppExport SEXP _NPMLEDIC_ep_DIC_vb(SEXP aind_LSEXP, SEXP aind_RSEXP, SEXP ctypeSEXP, SEXP alpha0SEXP, SEXP maxitSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat >::type aind_L(aind_LSEXP);
    Rcpp::traits::input_parameter< const arma::umat >::type aind_R(aind_RSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ctype(ctypeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(ep_DIC_vb(aind_L, aind_R, ctype, alpha0, maxit, tol));
    return rcpp_result_gen;
END_RCPP
}
// ep_DIC_gibbs
List ep_DIC_gibbs(const arma::umat aind_L, const arma::umat aind_R, const arma::uvec& ctype, const arma::vec& alpha0, const int& iter);
RcppExport SEXP _NPMLEDIC_ep_DIC_gibbs(SEXP aind_LSEXP, SEXP aind_RSEXP, SEXP ctypeSEXP, SEXP alpha0SEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat >::type aind_L(aind_LSEXP);
    Rcpp::traits::input_parameter< const arma::umat >::type aind_R(aind_RSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ctype(ctypeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(ep_DIC_gibbs(aind_L, aind_R, ctype, alpha0, iter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NPMLEDIC_ep_DIC_em", (DL_FUNC) &_NPMLEDIC_ep_DIC_em, 6},
    {"_NPMLEDIC_ep_DIC_vb", (DL_FUNC) &_NPMLEDIC_ep_DIC_vb, 6},
    {"_NPMLEDIC_ep_DIC_gibbs", (DL_FUNC) &_NPMLEDIC_ep_DIC_gibbs, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_NPMLEDIC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
