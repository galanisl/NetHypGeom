// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// hypermap
DataFrame hypermap(DataFrame net, double gma, double T, int k_speedup, double m_in, double L_in, double window, DataFrame theta);
RcppExport SEXP NetHypGeom_hypermap(SEXP netSEXP, SEXP gmaSEXP, SEXP TSEXP, SEXP k_speedupSEXP, SEXP m_inSEXP, SEXP L_inSEXP, SEXP windowSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type net(netSEXP);
    Rcpp::traits::input_parameter< double >::type gma(gmaSEXP);
    Rcpp::traits::input_parameter< double >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type k_speedup(k_speedupSEXP);
    Rcpp::traits::input_parameter< double >::type m_in(m_inSEXP);
    Rcpp::traits::input_parameter< double >::type L_in(L_inSEXP);
    Rcpp::traits::input_parameter< double >::type window(windowSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(hypermap(net, gma, T, k_speedup, m_in, L_in, window, theta));
    return rcpp_result_gen;
END_RCPP
}
