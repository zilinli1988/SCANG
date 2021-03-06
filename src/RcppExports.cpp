// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// matsp
List matsp(arma::mat G);
RcppExport SEXP _SCANG_matsp(SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(matsp(G));
    return rcpp_result_gen;
END_RCPP
}
// Liumod
double Liumod(arma::mat Cov, double q);
RcppExport SEXP _SCANG_Liumod(SEXP CovSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Cov(CovSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(Liumod(Cov, q));
    return rcpp_result_gen;
END_RCPP
}
// SCANG_O_Search
List SCANG_O_Search(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, arma::vec residuals, const double threshold_o, const double threshold_s, const double threshold_b, const int Lmax, const int Lmin, int steplength, arma::mat weights_B, arma::mat weights_S, const int begid, const double filter, const double f);
RcppExport SEXP _SCANG_SCANG_O_Search(SEXP GSEXP, SEXP XSEXP, SEXP workingSEXP, SEXP sigmaSEXP, SEXP famSEXP, SEXP residualsSEXP, SEXP threshold_oSEXP, SEXP threshold_sSEXP, SEXP threshold_bSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP steplengthSEXP, SEXP weights_BSEXP, SEXP weights_SSEXP, SEXP begidSEXP, SEXP filterSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type working(workingSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type fam(famSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type residuals(residualsSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_o(threshold_oSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_s(threshold_sSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_b(threshold_bSEXP);
    Rcpp::traits::input_parameter< const int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< const int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< int >::type steplength(steplengthSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_B(weights_BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_S(weights_SSEXP);
    Rcpp::traits::input_parameter< const int >::type begid(begidSEXP);
    Rcpp::traits::input_parameter< const double >::type filter(filterSEXP);
    Rcpp::traits::input_parameter< const double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(SCANG_O_Search(G, X, working, sigma, fam, residuals, threshold_o, threshold_s, threshold_b, Lmax, Lmin, steplength, weights_B, weights_S, begid, filter, f));
    return rcpp_result_gen;
END_RCPP
}
// SCANG_O_Search_Relatedness
List SCANG_O_Search_Relatedness(arma::sp_mat G, arma::mat P, arma::vec residuals, const double threshold_o, const double threshold_s, const double threshold_b, const int Lmax, const int Lmin, int steplength, arma::mat weights_B, arma::mat weights_S, const int begid, const double filter, const double f);
RcppExport SEXP _SCANG_SCANG_O_Search_Relatedness(SEXP GSEXP, SEXP PSEXP, SEXP residualsSEXP, SEXP threshold_oSEXP, SEXP threshold_sSEXP, SEXP threshold_bSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP steplengthSEXP, SEXP weights_BSEXP, SEXP weights_SSEXP, SEXP begidSEXP, SEXP filterSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type residuals(residualsSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_o(threshold_oSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_s(threshold_sSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_b(threshold_bSEXP);
    Rcpp::traits::input_parameter< const int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< const int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< int >::type steplength(steplengthSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_B(weights_BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_S(weights_SSEXP);
    Rcpp::traits::input_parameter< const int >::type begid(begidSEXP);
    Rcpp::traits::input_parameter< const double >::type filter(filterSEXP);
    Rcpp::traits::input_parameter< const double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(SCANG_O_Search_Relatedness(G, P, residuals, threshold_o, threshold_s, threshold_b, Lmax, Lmin, steplength, weights_B, weights_S, begid, filter, f));
    return rcpp_result_gen;
END_RCPP
}
// SCANG_O_Search_Relatedness_sp
List SCANG_O_Search_Relatedness_sp(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::vec residuals, const double threshold_o, const double threshold_s, const double threshold_b, const int Lmax, const int Lmin, int steplength, arma::mat weights_B, arma::mat weights_S, const int begid, const double filter, const double f);
RcppExport SEXP _SCANG_SCANG_O_Search_Relatedness_sp(SEXP GSEXP, SEXP Sigma_iSEXP, SEXP Sigma_iXSEXP, SEXP covSEXP, SEXP residualsSEXP, SEXP threshold_oSEXP, SEXP threshold_sSEXP, SEXP threshold_bSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP steplengthSEXP, SEXP weights_BSEXP, SEXP weights_SSEXP, SEXP begidSEXP, SEXP filterSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Sigma_i(Sigma_iSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_iX(Sigma_iXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov(covSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type residuals(residualsSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_o(threshold_oSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_s(threshold_sSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_b(threshold_bSEXP);
    Rcpp::traits::input_parameter< const int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< const int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< int >::type steplength(steplengthSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_B(weights_BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_S(weights_SSEXP);
    Rcpp::traits::input_parameter< const int >::type begid(begidSEXP);
    Rcpp::traits::input_parameter< const double >::type filter(filterSEXP);
    Rcpp::traits::input_parameter< const double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(SCANG_O_Search_Relatedness_sp(G, Sigma_i, Sigma_iX, cov, residuals, threshold_o, threshold_s, threshold_b, Lmax, Lmin, steplength, weights_B, weights_S, begid, filter, f));
    return rcpp_result_gen;
END_RCPP
}
// SCANG_O_Thres
arma::mat SCANG_O_Thres(arma::sp_mat G, arma::mat X, arma::vec working, double sigma, int fam, int times, int Lmax, int Lmin, int steplength, arma::mat weights_B, arma::mat weights_S, double filter);
RcppExport SEXP _SCANG_SCANG_O_Thres(SEXP GSEXP, SEXP XSEXP, SEXP workingSEXP, SEXP sigmaSEXP, SEXP famSEXP, SEXP timesSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP steplengthSEXP, SEXP weights_BSEXP, SEXP weights_SSEXP, SEXP filterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type working(workingSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type fam(famSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    Rcpp::traits::input_parameter< int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< int >::type steplength(steplengthSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_B(weights_BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_S(weights_SSEXP);
    Rcpp::traits::input_parameter< double >::type filter(filterSEXP);
    rcpp_result_gen = Rcpp::wrap(SCANG_O_Thres(G, X, working, sigma, fam, times, Lmax, Lmin, steplength, weights_B, weights_S, filter));
    return rcpp_result_gen;
END_RCPP
}
// SCANG_O_Thres_Relatedness
arma::mat SCANG_O_Thres_Relatedness(arma::sp_mat G, arma::mat P, arma::mat pesudo_residuals, int times, int Lmax, int Lmin, int steplength, arma::mat weights_B, arma::mat weights_S, double filter);
RcppExport SEXP _SCANG_SCANG_O_Thres_Relatedness(SEXP GSEXP, SEXP PSEXP, SEXP pesudo_residualsSEXP, SEXP timesSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP steplengthSEXP, SEXP weights_BSEXP, SEXP weights_SSEXP, SEXP filterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pesudo_residuals(pesudo_residualsSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    Rcpp::traits::input_parameter< int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< int >::type steplength(steplengthSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_B(weights_BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_S(weights_SSEXP);
    Rcpp::traits::input_parameter< double >::type filter(filterSEXP);
    rcpp_result_gen = Rcpp::wrap(SCANG_O_Thres_Relatedness(G, P, pesudo_residuals, times, Lmax, Lmin, steplength, weights_B, weights_S, filter));
    return rcpp_result_gen;
END_RCPP
}
// SCANG_O_Thres_Relatedness_sp
arma::mat SCANG_O_Thres_Relatedness_sp(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::mat pesudo_residuals, int times, int Lmax, int Lmin, int steplength, arma::mat weights_B, arma::mat weights_S, double filter);
RcppExport SEXP _SCANG_SCANG_O_Thres_Relatedness_sp(SEXP GSEXP, SEXP Sigma_iSEXP, SEXP Sigma_iXSEXP, SEXP covSEXP, SEXP pesudo_residualsSEXP, SEXP timesSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP steplengthSEXP, SEXP weights_BSEXP, SEXP weights_SSEXP, SEXP filterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Sigma_i(Sigma_iSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_iX(Sigma_iXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov(covSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pesudo_residuals(pesudo_residualsSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    Rcpp::traits::input_parameter< int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< int >::type steplength(steplengthSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_B(weights_BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_S(weights_SSEXP);
    Rcpp::traits::input_parameter< double >::type filter(filterSEXP);
    rcpp_result_gen = Rcpp::wrap(SCANG_O_Thres_Relatedness_sp(G, Sigma_i, Sigma_iX, cov, pesudo_residuals, times, Lmax, Lmin, steplength, weights_B, weights_S, filter));
    return rcpp_result_gen;
END_RCPP
}
// maxO
arma::mat maxO(int p, int Lmax, int Lmin, int steplength, arma::mat x, arma::mat weights_B, arma::mat weights_S, arma::mat Cov, double filter, int times);
RcppExport SEXP _SCANG_maxO(SEXP pSEXP, SEXP LmaxSEXP, SEXP LminSEXP, SEXP steplengthSEXP, SEXP xSEXP, SEXP weights_BSEXP, SEXP weights_SSEXP, SEXP CovSEXP, SEXP filterSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type Lmax(LmaxSEXP);
    Rcpp::traits::input_parameter< int >::type Lmin(LminSEXP);
    Rcpp::traits::input_parameter< int >::type steplength(steplengthSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_B(weights_BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weights_S(weights_SSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cov(CovSEXP);
    Rcpp::traits::input_parameter< double >::type filter(filterSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(maxO(p, Lmax, Lmin, steplength, x, weights_B, weights_S, Cov, filter, times));
    return rcpp_result_gen;
END_RCPP
}
// regionfilter
arma::mat regionfilter(arma::mat candidate, const double f);
RcppExport SEXP _SCANG_regionfilter(SEXP candidateSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type candidate(candidateSEXP);
    Rcpp::traits::input_parameter< const double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(regionfilter(candidate, f));
    return rcpp_result_gen;
END_RCPP
}
// svd_c
List svd_c(arma::mat X);
RcppExport SEXP _SCANG_svd_c(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(svd_c(X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SCANG_matsp", (DL_FUNC) &_SCANG_matsp, 1},
    {"_SCANG_Liumod", (DL_FUNC) &_SCANG_Liumod, 2},
    {"_SCANG_SCANG_O_Search", (DL_FUNC) &_SCANG_SCANG_O_Search, 17},
    {"_SCANG_SCANG_O_Search_Relatedness", (DL_FUNC) &_SCANG_SCANG_O_Search_Relatedness, 14},
    {"_SCANG_SCANG_O_Search_Relatedness_sp", (DL_FUNC) &_SCANG_SCANG_O_Search_Relatedness_sp, 16},
    {"_SCANG_SCANG_O_Thres", (DL_FUNC) &_SCANG_SCANG_O_Thres, 12},
    {"_SCANG_SCANG_O_Thres_Relatedness", (DL_FUNC) &_SCANG_SCANG_O_Thres_Relatedness, 10},
    {"_SCANG_SCANG_O_Thres_Relatedness_sp", (DL_FUNC) &_SCANG_SCANG_O_Thres_Relatedness_sp, 12},
    {"_SCANG_maxO", (DL_FUNC) &_SCANG_maxO, 10},
    {"_SCANG_regionfilter", (DL_FUNC) &_SCANG_regionfilter, 2},
    {"_SCANG_svd_c", (DL_FUNC) &_SCANG_svd_c, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_SCANG(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
