// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// disj
arma::umat disj(arma::uvec strata);
RcppExport SEXP _fastcube_disj(SEXP strataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type strata(strataSEXP);
    rcpp_result_gen = Rcpp::wrap(disj(strata));
    return rcpp_result_gen;
END_RCPP
}
// ncat
arma::rowvec ncat(arma::umat Xcat);
RcppExport SEXP _fastcube_ncat(SEXP XcatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type Xcat(XcatSEXP);
    rcpp_result_gen = Rcpp::wrap(ncat(Xcat));
    return rcpp_result_gen;
END_RCPP
}
// disjMatrix
arma::umat disjMatrix(arma::umat strata);
RcppExport SEXP _fastcube_disjMatrix(SEXP strataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type strata(strataSEXP);
    rcpp_result_gen = Rcpp::wrap(disjMatrix(strata));
    return rcpp_result_gen;
END_RCPP
}
// findBarma
arma::mat findBarma(arma::mat X, arma::umat Xcat);
RcppExport SEXP _fastcube_findBarma(SEXP XSEXP, SEXP XcatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type Xcat(XcatSEXP);
    rcpp_result_gen = Rcpp::wrap(findBarma(X, Xcat));
    return rcpp_result_gen;
END_RCPP
}
// isEye
bool isEye(arma::mat& M);
RcppExport SEXP _fastcube_isEye(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(isEye(M));
    return rcpp_result_gen;
END_RCPP
}
// rrefArma
void rrefArma(arma::mat& M);
RcppExport SEXP _fastcube_rrefArma(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    rrefArma(M);
    return R_NilValue;
END_RCPP
}
// osffphase
arma::vec osffphase(arma::vec prob, arma::mat Bm);
RcppExport SEXP _fastcube_osffphase(SEXP probSEXP, SEXP BmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type prob(probSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Bm(BmSEXP);
    rcpp_result_gen = Rcpp::wrap(osffphase(prob, Bm));
    return rcpp_result_gen;
END_RCPP
}
// onestep
arma::vec onestep(arma::mat B, arma::vec pik, double EPS);
RcppExport SEXP _fastcube_onestep(SEXP BSEXP, SEXP pikSEXP, SEXP EPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< double >::type EPS(EPSSEXP);
    rcpp_result_gen = Rcpp::wrap(onestep(B, pik, EPS));
    return rcpp_result_gen;
END_RCPP
}
// flightphase_arma
arma::vec flightphase_arma(arma::mat X, arma::vec pik, bool redux);
RcppExport SEXP _fastcube_flightphase_arma(SEXP XSEXP, SEXP pikSEXP, SEXP reduxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< bool >::type redux(reduxSEXP);
    rcpp_result_gen = Rcpp::wrap(flightphase_arma(X, pik, redux));
    return rcpp_result_gen;
END_RCPP
}
// reduxArma
Rcpp::List reduxArma(arma::mat B);
RcppExport SEXP _fastcube_reduxArma(SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(reduxArma(B));
    return rcpp_result_gen;
END_RCPP
}
// systematicDesign
Rcpp::List systematicDesign(arma::vec pik);
RcppExport SEXP _fastcube_systematicDesign(SEXP pikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    rcpp_result_gen = Rcpp::wrap(systematicDesign(pik));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastcube_disj", (DL_FUNC) &_fastcube_disj, 1},
    {"_fastcube_ncat", (DL_FUNC) &_fastcube_ncat, 1},
    {"_fastcube_disjMatrix", (DL_FUNC) &_fastcube_disjMatrix, 1},
    {"_fastcube_findBarma", (DL_FUNC) &_fastcube_findBarma, 2},
    {"_fastcube_isEye", (DL_FUNC) &_fastcube_isEye, 1},
    {"_fastcube_rrefArma", (DL_FUNC) &_fastcube_rrefArma, 1},
    {"_fastcube_osffphase", (DL_FUNC) &_fastcube_osffphase, 2},
    {"_fastcube_onestep", (DL_FUNC) &_fastcube_onestep, 3},
    {"_fastcube_flightphase_arma", (DL_FUNC) &_fastcube_flightphase_arma, 3},
    {"_fastcube_reduxArma", (DL_FUNC) &_fastcube_reduxArma, 1},
    {"_fastcube_systematicDesign", (DL_FUNC) &_fastcube_systematicDesign, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastcube(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
