// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dispen_map
NumericVector dispen_map(NumericVector x, int ma, int nc, double mu, double sigma);
RcppExport SEXP _MSEntropy_dispen_map(SEXP xSEXP, SEXP maSEXP, SEXP ncSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ma(maSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(dispen_map(x, ma, nc, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// dispen_npdf
NumericVector dispen_npdf(NumericVector z, int nc, int m, int tau);
RcppExport SEXP _MSEntropy_dispen_npdf(SEXP zSEXP, SEXP ncSEXP, SEXP mSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(dispen_npdf(z, nc, m, tau));
    return rcpp_result_gen;
END_RCPP
}
// fdispen_npdf
NumericVector fdispen_npdf(NumericVector z, int nc, int m, int tau);
RcppExport SEXP _MSEntropy_fdispen_npdf(SEXP zSEXP, SEXP ncSEXP, SEXP mSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(fdispen_npdf(z, nc, m, tau));
    return rcpp_result_gen;
END_RCPP
}
// SampEn
double SampEn(NumericVector x, int m, double r, double sd);
RcppExport SEXP _MSEntropy_SampEn(SEXP xSEXP, SEXP mSEXP, SEXP rSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(SampEn(x, m, r, sd));
    return rcpp_result_gen;
END_RCPP
}
// CoarseGraining
NumericVector CoarseGraining(NumericVector x, int scale);
RcppExport SEXP _MSEntropy_CoarseGraining(SEXP xSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(CoarseGraining(x, scale));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MSEntropy_dispen_map", (DL_FUNC) &_MSEntropy_dispen_map, 5},
    {"_MSEntropy_dispen_npdf", (DL_FUNC) &_MSEntropy_dispen_npdf, 4},
    {"_MSEntropy_fdispen_npdf", (DL_FUNC) &_MSEntropy_fdispen_npdf, 4},
    {"_MSEntropy_SampEn", (DL_FUNC) &_MSEntropy_SampEn, 4},
    {"_MSEntropy_CoarseGraining", (DL_FUNC) &_MSEntropy_CoarseGraining, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MSEntropy(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
