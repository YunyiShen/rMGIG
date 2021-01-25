// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// CARE_ArimotoPotter_cpp
void CARE_ArimotoPotter_cpp(arma::mat& X, const arma::mat& A, const arma::mat& R, const arma::mat& Q);
RcppExport SEXP _rMGIG_CARE_ArimotoPotter_cpp(SEXP XSEXP, SEXP ASEXP, SEXP RSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    CARE_ArimotoPotter_cpp(X, A, R, Q);
    return R_NilValue;
END_RCPP
}
// rMGIG_cpp
arma::mat rMGIG_cpp(int n, const double& nu, const arma::mat& phi, const arma::mat& psi, const double& df, int maxit);
RcppExport SEXP _rMGIG_rMGIG_cpp(SEXP nSEXP, SEXP nuSEXP, SEXP phiSEXP, SEXP psiSEXP, SEXP dfSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(rMGIG_cpp(n, nu, phi, psi, df, maxit));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rMGIG_CARE_ArimotoPotter_cpp", (DL_FUNC) &_rMGIG_CARE_ArimotoPotter_cpp, 4},
    {"_rMGIG_rMGIG_cpp", (DL_FUNC) &_rMGIG_rMGIG_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_rMGIG(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
