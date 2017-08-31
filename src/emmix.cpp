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


Rcpp::NumericVector export_vec(arma::vec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}

Rcpp::NumericVector export_uvec(arma::uvec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
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

arma::mat start_kmeans(arma::vec dat, int g){
  arma::mat params(g,4, arma::fill::zeros);
  
  //import kmeans
  
  return(params);
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


double BIC_calc(double LL,  int k, int n){
  double BIC = log(n)* k - 2*LL ;
  return(BIC);
}

double AIC_calc(double LL,  int k){
  double AIC = 2*k - 2*LL;
  return(AIC);
}



// [[Rcpp::export]]
List emmix_t(arma::vec dat, int g=1, int random_starts=4, int max_it=100, double tol = 0.0001){
  int n = dat.size();
  int k =0; int n_fixed =0;
  double  LL =0;
  arma::mat params(g,4, arma::fill::zeros);
  
  arma::vec pi = arma::zeros(g);
  arma::vec mu = arma::zeros(g);
  arma::vec nu = arma::zeros(g);
  arma::vec sigma = arma::zeros(g);
 
  arma::uvec clusters = arma::zeros<arma::uvec>(n);
  arma::uvec cluster_sizes = arma::zeros<arma::uvec>(g);
  arma::vec lik = arma::zeros(max_it);
  
  arma::mat tau = arma::zeros(g,n);
  arma::mat you = arma::zeros(g,n);
  arma::mat temp2(g,4, arma::fill::zeros);
  
  double Q1 = 0.0;
  arma::vec Q2 =  arma::zeros(g,n);
  arma::vec Q3 =  arma::zeros(g,n);

  params = start_kmeans(dat, g);
  
  //primary loop
  
  for(int j =0; j<max_it;j++){
    arma::mat temp1 = estep(dat,params);
    tau = temp1.rows(0,g-1);
    you = temp1.rows(g,2*g-1);
    temp2 = mstep(dat, tau, you, params);
    params = temp2; 
    //log liklihood
    Q1 = sum(sum( tau.each_col() % arma::log(params.col(0)) ));
   
    mu=params.col(1);
    nu=params.col(2);
    sigma = params.col(3);
    
    for(int i=0; i <g; i++){
      Q2.row(i) = -log ( gamma(0.5*(nu.at(i)))   ) + 0.5*nu.at(i) *log(0.5*nu.at(i))+ 0.5*nu.at(i)*sum(log(you.row(i))-you.row(i),1) + digamma((nu.at(i)+1)/2) - log ((nu.at(i)+1)/2);
      Q3.row(i) = -0.5*1*log(2*arma::datum::pi) - 0.5*log(abs(sigma.at(i))) + 0.5*1*log(you.row(i)) - 0.5*you.row(i)*(1/sigma.at(i))*((dat-mu.at(i)).t())*(dat-mu.at(i));
    }
    
    lik.at(j) = Q1 + sum(sum(tau % Q2,0)) + sum(sum(tau % Q3,0));
    LL = lik.at(j);
    
    
  }
  
  //cluster allocation
  for(int i =0; i<n;i++){
    clusters.at(i) = tau.col(i).index_max();
  }
  
  //cluster size
  for(int i =0; i<g;i++){
    arma::uvec tmp1 = find(clusters == i);
    cluster_sizes.at(i) = tmp1.n_elem;
  }
  //min cluster size
  int c_min = min(cluster_sizes);
  
  
  //n_params
  k = (g-1)+ g + g + g - n_fixed;//pi + mu + sigma + nu - fixed params
  
  return  List::create(
                      Named("mu")= export_vec(params.col(1)),
                      Named("pi") = export_vec(params.col(0)),
                      Named("nu") = export_vec(params.col(2)),
                      Named("sigma") = export_vec(params.col(3)),
                      Named("LL") = LL,
                      Named("BIC") = BIC_calc(LL, k, n),
                      Named("AIC") = AIC_calc(LL, k),
                      Named("tau") = tau,
                      Named("lik") = lik,
                      Named("k") = k, //number of parameters
                      Named("Clusters") = export_uvec(clusters),
                      Named("c_min") = c_min //smallest cluster
                      );

}


List each_gene(arma::vec dat, int random_starts=4, double ll_thresh = 8, int min_clust_size = 8, double tol = 0.0001){
  
  List g1 = emmix_t(dat, 1);
  List g2 = emmix_t(dat, 2);
  List best_g = g1;
  
  double lambda = 0; // ll_ratio(List g1, List g2)
  
  if(-2*log(lambda) > ll_thresh){
    if(as<int>(g2.attr("c_min")) > min_clust_size){
      best_g = g2;
    }else{
      List g3 = emmix_t(dat, 3, 100);
      lambda = 0; // ll_ratio(List g2, List g3)
      if(-2*log(lambda) > ll_thresh){
        if(as<int>(g3.attr("c_min")) > min_clust_size){
          best_g = g3;
        }
      }
    }
  }
  
  return(best_g);
}

// [[Rcpp::export]]
List emmix_gene(arma::mat bigdat, int random_starts=4, double ll_thresh = 8, int min_clust_size = 8, double tol = 0.0001){
 
 int n = bigdat.n_rows;
 List tmp = List::create();
 
 for(int i; i<n;i++){
   tmp = each_gene(bigdat.row(i), random_starts, ll_thresh, min_clust_size, tol);
 }
 
 



  return(tmp);
}












