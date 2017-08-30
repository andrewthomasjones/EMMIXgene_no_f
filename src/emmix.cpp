//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//'@importFrom Rcpp sourceCpp

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <algorithm>
#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/students_t.hpp>

//using namespace arma;
using namespace Rcpp;
using namespace boost::math;

typedef students_t_distribution<double> students_t;
typedef std::vector<double> stdvec;

double df_eq_est (double v, double v_sum);
//function class with operator for root finder for df of t-dist
class df_eq_func{
public:
  df_eq_func(double v) :  v_sum(v){};

  double getSumV();

  double operator()(double v) {

    double temp = -digamma(0.5*v) + std::log(0.5*v) + v_sum + digamma((v+1.)/2.) - std::log((v+1.)/2.)+1.;
    //std::cout << "v: " << v << " v_sum: " << v_sum << " temp: " << temp << std::endl;
    return temp;
    //return df_eq_est(v,v_sum);

  }

private:
  double v_sum;
};

double df_eq_func::getSumV()
{
  return v_sum;
}


// [[Rcpp::export]]
arma::mat estep(const arma::vec& dat, const arma::mat& params){
  int n = dat.size();
  int g = params.n_rows;
  int g2 =2*g;
  
  arma::mat new_params(g2, n, arma::fill::zeros);
  arma::mat tau(g,n, arma::fill::zeros);
  arma::mat tau_sum(1,n, arma::fill::zeros);
  arma::mat you(g,n, arma::fill::zeros);

  arma::vec pi = params.col(0);
  arma::vec nu = params.col(2);
  arma::vec mu = params.col(1);
  arma::vec sigma = params.col(3);


  for(int i=0;i<g;i++){
    students_t dist(nu(i));
    stdvec dens(n);
    std::transform(dat.begin(), dat.end(), dens.begin(), [dist](double dat) { return pdf(dist,dat); });
    arma::vec dens2 = arma::conv_to<arma::vec>::from(dens);
    tau.row(i) = pi.at(i)*dens2.t();
    tau_sum.row(0) += pi.at(i)*dens2.t();
    you.row(i) = ((nu.at(i)+1.0)/(nu.at(i) + (dat-mu.at(i))%(dat-mu.at(i))*(1/sigma.at(i)) )).t();

  }
  
  for(int i=0;i<g;i++){
    tau.row(i) /= tau_sum.row(0);
  }

  new_params.rows(0,g-1)=tau;
  new_params.rows(g,g2-1)=you;

  return(new_params);
}




// [[Rcpp::export]]
arma::mat mstep(arma::vec dat, arma::mat tau, arma::mat you, arma::mat params){
  int n = dat.size();
  int g = tau.n_rows;

  arma::vec mu2 = arma::zeros<arma::vec>(g);
  arma::vec nu2 = arma::zeros<arma::vec>(g);
  arma::vec sigma2 = arma::zeros<arma::vec>(g);
  arma::vec s_tau = arma::zeros<arma::vec>(g);
  arma::vec pi2 = arma::zeros<arma::vec>(g);
  
  s_tau = sum(tau,1);
  pi2 = s_tau/n;

  double a = 0.001;
  double b = 100;
  boost::uintmax_t df_max_iter=500; // boost solver params, could be user params but relatively unimportant
  tools::eps_tolerance<double> tol(30);

  for(int i=0;i<g;i++){
    mu2(i) = (sum(dat.t()%you.row(i)%tau.row(i)) / sum(you.row(i)%tau.row(i)));
    sigma2(i) = (sum(you.row(i)%tau.row(i)%((dat.t()-mu2(i))%(dat.t()-mu2(i)))  ) / s_tau(i));

    double v_sum = (1/s_tau(i))*sum(tau.row(i)%(arma::log(you.row(i))-you.row(i)));
    df_eq_func rootFun = df_eq_func(v_sum);
    std::pair<double, double>  r1= tools::bisect(rootFun, a, b, tol, df_max_iter);
    nu2(i) = (r1.first + r1.second)/2.0;

  }

  arma::mat new_params(g,4, arma::fill::zeros);
  new_params.col(0) = pi2;
  new_params.col(1) = mu2;
  new_params.col(2) = nu2;
  new_params.col(3) = sigma2;

  return(new_params);
}


// [[Rcpp::export]]
List emmix_t(arma::vec dat, arma::mat params, int g=1, int max_it=100){
  int n = dat.size();

  arma::vec pi = arma::zeros(g);
  arma::vec mu = arma::zeros(g);
  arma::vec nu = arma::zeros(g);
  arma::vec sigma = arma::zeros(g);
 
  arma::uvec clusters; clusters.fill(0);
  arma::vec lik = arma::zeros(max_it);
  
  arma::mat tau = arma::zeros(g,n);
  arma::mat you = arma::zeros(g,n);
  arma::mat temp2(g,4, arma::fill::zeros);
  
  double Q1 =0.0;

  arma::vec Q2 = arma::zeros(n);
  arma::vec Q3 = arma::zeros(n);


  pi =  params.col(0);
  mu =  params.col(1);
  nu =  params.col(2);
  sigma =  params.col(3);

  
  //primary loop
  
  for(int k =0; k<max_it;k++){
    arma::mat temp1 = estep(dat,params);
    tau = temp1.rows(0,g-1);
    you = temp1.rows(g,2*g-1);
    temp2 = mstep(dat, tau, you, params);
  }
  
  params = temp2;  
  
 
  
  
  return Rcpp::List::create(Rcpp::Named("mu")= params.col(1));
                      //   
                      // Named("pi") = params.col(0),
                      // Named("nu") = params.col(2),
                      // Named("sigma") = params.col(3),
                      // Named("LL") = 0,
                      // Named("BIC") = 0,  
                      // Named("AIC") = 0,
                      // Named("tau") = tau,
                      // Named("lik") = lik,
                      // Named("Clusters") = clusters
                      // );

}





















