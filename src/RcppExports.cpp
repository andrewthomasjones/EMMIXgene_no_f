// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// call_emmix_sel
arma:: vec call_emmix_sel(arma:: vec t, int row, int col, int g);
RcppExport SEXP _EMMIXgene_call_emmix_sel(SEXP tSEXP, SEXP rowSEXP, SEXP colSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma:: vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type row(rowSEXP);
    Rcpp::traits::input_parameter< int >::type col(colSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(call_emmix_sel(t, row, col, g));
    return rcpp_result_gen;
END_RCPP
}
// select_genes
int select_genes(arma::mat& data, int row, int col, int g, int k, int r, double b1, int b2);
RcppExport SEXP _EMMIXgene_select_genes(SEXP dataSEXP, SEXP rowSEXP, SEXP colSEXP, SEXP gSEXP, SEXP kSEXP, SEXP rSEXP, SEXP b1SEXP, SEXP b2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type row(rowSEXP);
    Rcpp::traits::input_parameter< int >::type col(colSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< int >::type b2(b2SEXP);
    rcpp_result_gen = Rcpp::wrap(select_genes(data, row, col, g, k, r, b1, b2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EMMIXgene_call_emmix_sel", (DL_FUNC) &_EMMIXgene_call_emmix_sel, 4},
    {"_EMMIXgene_select_genes", (DL_FUNC) &_EMMIXgene_select_genes, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_EMMIXgene(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
