// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mahalanobis
double mahalanobis(double y, double mu, double sigma);
RcppExport SEXP _EMMIXgene_mahalanobis(SEXP ySEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mahalanobis(y, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// t_dist
double t_dist(double y, double mu, double sigma, double nu, int p);
RcppExport SEXP _EMMIXgene_t_dist(SEXP ySEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nuSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(t_dist(y, mu, sigma, nu, p));
    return rcpp_result_gen;
END_RCPP
}
// estep
arma::mat estep(const arma::vec& dat, const arma::mat& params);
RcppExport SEXP _EMMIXgene_estep(SEXP datSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type dat(datSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(estep(dat, params));
    return rcpp_result_gen;
END_RCPP
}
// mstep
arma::mat mstep(arma::vec dat, arma::mat tau, arma::mat you, arma::mat params);
RcppExport SEXP _EMMIXgene_mstep(SEXP datSEXP, SEXP tauSEXP, SEXP youSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type you(youSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(mstep(dat, tau, you, params));
    return rcpp_result_gen;
END_RCPP
}
// emmix_t
List emmix_t(arma::vec dat, int g, int random_starts, int max_it, double tol, std::string start_method);
RcppExport SEXP _EMMIXgene_emmix_t(SEXP datSEXP, SEXP gSEXP, SEXP random_startsSEXP, SEXP max_itSEXP, SEXP tolSEXP, SEXP start_methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type random_starts(random_startsSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< std::string >::type start_method(start_methodSEXP);
    rcpp_result_gen = Rcpp::wrap(emmix_t(dat, g, random_starts, max_it, tol, start_method));
    return rcpp_result_gen;
END_RCPP
}
// each_gene
List each_gene(arma::vec dat, int random_starts, int max_it, double ll_thresh, int min_clust_size, double tol, std::string start_method);
RcppExport SEXP _EMMIXgene_each_gene(SEXP datSEXP, SEXP random_startsSEXP, SEXP max_itSEXP, SEXP ll_threshSEXP, SEXP min_clust_sizeSEXP, SEXP tolSEXP, SEXP start_methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type random_starts(random_startsSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type ll_thresh(ll_threshSEXP);
    Rcpp::traits::input_parameter< int >::type min_clust_size(min_clust_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< std::string >::type start_method(start_methodSEXP);
    rcpp_result_gen = Rcpp::wrap(each_gene(dat, random_starts, max_it, ll_thresh, min_clust_size, tol, start_method));
    return rcpp_result_gen;
END_RCPP
}
// emmix_gene
List emmix_gene(arma::mat& bigdat, int random_starts, int max_it, double ll_thresh, int min_clust_size, double tol, std::string start_method);
RcppExport SEXP _EMMIXgene_emmix_gene(SEXP bigdatSEXP, SEXP random_startsSEXP, SEXP max_itSEXP, SEXP ll_threshSEXP, SEXP min_clust_sizeSEXP, SEXP tolSEXP, SEXP start_methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type bigdat(bigdatSEXP);
    Rcpp::traits::input_parameter< int >::type random_starts(random_startsSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type ll_thresh(ll_threshSEXP);
    Rcpp::traits::input_parameter< int >::type min_clust_size(min_clust_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< std::string >::type start_method(start_methodSEXP);
    rcpp_result_gen = Rcpp::wrap(emmix_gene(bigdat, random_starts, max_it, ll_thresh, min_clust_size, tol, start_method));
    return rcpp_result_gen;
END_RCPP
}
// tkmeans
arma::mat tkmeans(arma::mat& M, int k, double alpha, arma::vec weights, int nstart, int iter, double tol, bool verbose);
RcppExport SEXP _EMMIXgene_tkmeans(SEXP MSEXP, SEXP kSEXP, SEXP alphaSEXP, SEXP weightsSEXP, SEXP nstartSEXP, SEXP iterSEXP, SEXP tolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type nstart(nstartSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(tkmeans(M, k, alpha, weights, nstart, iter, tol, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EMMIXgene_mahalanobis", (DL_FUNC) &_EMMIXgene_mahalanobis, 3},
    {"_EMMIXgene_t_dist", (DL_FUNC) &_EMMIXgene_t_dist, 5},
    {"_EMMIXgene_estep", (DL_FUNC) &_EMMIXgene_estep, 2},
    {"_EMMIXgene_mstep", (DL_FUNC) &_EMMIXgene_mstep, 4},
    {"_EMMIXgene_emmix_t", (DL_FUNC) &_EMMIXgene_emmix_t, 6},
    {"_EMMIXgene_each_gene", (DL_FUNC) &_EMMIXgene_each_gene, 7},
    {"_EMMIXgene_emmix_gene", (DL_FUNC) &_EMMIXgene_emmix_gene, 7},
    {"_EMMIXgene_tkmeans", (DL_FUNC) &_EMMIXgene_tkmeans, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_EMMIXgene(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
