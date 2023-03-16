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
List ep_DIC_em(const arma::vec& EL, const arma::vec& ER, const arma::vec& SL, const arma::vec& SR, const arma::uvec& ctype, const arma::vec& breaks, const double& alpha0, const int& iter);
RcppExport SEXP _NPMLEDIC_ep_DIC_em(SEXP ELSEXP, SEXP ERSEXP, SEXP SLSEXP, SEXP SRSEXP, SEXP ctypeSEXP, SEXP breaksSEXP, SEXP alpha0SEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type EL(ELSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ER(ERSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type SL(SLSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type SR(SRSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ctype(ctypeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type breaks(breaksSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(ep_DIC_em(EL, ER, SL, SR, ctype, breaks, alpha0, iter));
    return rcpp_result_gen;
END_RCPP
}
// ep_DIC_vb
List ep_DIC_vb(const arma::vec& EL, const arma::vec& ER, const arma::vec& SL, const arma::vec& SR, const arma::uvec& ctype, const arma::vec& breaks, const double& alpha0, const int& iter);
RcppExport SEXP _NPMLEDIC_ep_DIC_vb(SEXP ELSEXP, SEXP ERSEXP, SEXP SLSEXP, SEXP SRSEXP, SEXP ctypeSEXP, SEXP breaksSEXP, SEXP alpha0SEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type EL(ELSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ER(ERSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type SL(SLSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type SR(SRSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ctype(ctypeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type breaks(breaksSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(ep_DIC_vb(EL, ER, SL, SR, ctype, breaks, alpha0, iter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NPMLEDIC_ep_DIC_em", (DL_FUNC) &_NPMLEDIC_ep_DIC_em, 8},
    {"_NPMLEDIC_ep_DIC_vb", (DL_FUNC) &_NPMLEDIC_ep_DIC_vb, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_NPMLEDIC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
