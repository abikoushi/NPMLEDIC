// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ep_ICRT_em
arma::vec ep_ICRT_em(const arma::vec& EL, const arma::vec& ER, const arma::vec& S, const double& tmax, const arma::vec& breaks, const int& iter);
RcppExport SEXP _NPMLEDIC_ep_ICRT_em(SEXP ELSEXP, SEXP ERSEXP, SEXP SSEXP, SEXP tmaxSEXP, SEXP breaksSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type EL(ELSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ER(ERSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double& >::type tmax(tmaxSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type breaks(breaksSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(ep_ICRT_em(EL, ER, S, tmax, breaks, iter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NPMLEDIC_ep_ICRT_em", (DL_FUNC) &_NPMLEDIC_ep_ICRT_em, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_NPMLEDIC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}