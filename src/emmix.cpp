//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//using namespace arma;
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//'@importFrom Rcpp sourceCpp
//'@useDynLib EMMIXgene

#include <vector>
#include <string>
#include <algorithm>

arma::vec estep(arma::vec dat, arma::vec params){



  return(params);
}

arma::vec mstep(arma::vec dat, arma::vec params){



  return(params);
}


// [[Rcpp::export]]
arma::vec emmix_t(arma::vec dat){
  int col = dat.size();


  return(dat);
}





















