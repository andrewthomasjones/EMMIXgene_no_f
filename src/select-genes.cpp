//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <emmix.h>

//using namespace arma;
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//'@importFrom Rcpp sourceCpp
//'@useDynLib EMMIXgene

#include <vector>
#include <string>
#include <algorithm>

#define LINE_LENGTH 5000	/* max length of each line in microarray data */
#define MAX_GENES 5000		/* max number of genes in microarray data */
#define MAX_TISSUES 200     /* max number of tissues in microarray data */
#define MAX_GROUPS 100		/* max number of groups to split into */
#define MAX_ENTRY_LENGTH 20	/* max length of any single entry in the array */
#define MAX_SELECTED_GENES 100



struct tab
{
  int nu;
  int count;
  int sm;
  double tl;
  std::vector<double> s;
};

bool tabcompare(const tab &a, const tab &b)
{
  return a.tl< b.tl;
};




//&tl, &sm, &er

// [[Rcpp::export]]
 arma::vec call_emmix_sel( arma::vec t, int row, int col, int g){
  double fl1, fl2;
  int group1[MAX_TISSUES], group2[MAX_TISSUES];
  int comp[MAX_TISSUES], x, y = 0, a1, a2, i;
   arma::vec s =  arma::zeros(col);
   arma::vec ret =  arma::zeros(3);
  double tl, sm, er;

  if (g == 1){
    //actually call emmix on t
  }
  else{
    //actually call emmix on t
  }

  //out1.%d"
  a1 = 0;//sscanf (s, " Final Log-Likelihood is %f", &fl1);
  //out2.%d
  a2 = 0;//sscanf (s, " Final Log-Likelihood is %f", &fl2);

  //final loglikes
  fl1=0;
  fl2=0;

  if (a1 && a2){
    tl = 2.0 * (fl2 - fl1);
  }
  else{
    ret(0) = tl;
    ret(1) = er;
    ret(2) = sm;
    return(ret);
  }


  //out2
  // group2[i] = some number per line

  er = 0;

  memset (comp, 0, sizeof (comp));
  for (i = 0; i < col; i++){
    comp[group2[i]]++;
  }

  sm = 9999;
  for (i = 1; i <= g + 1; i++){
    if (comp[i] < sm){
      sm = comp[i];
    }
  }

  /* 06/04/2002 CGT GJM
   ** If minimum size is 1, the set -2 log lambda t = -100
   */
  if (sm == 1)
  {
    tl = -100;
    er = -100;
  }

  //return stuff
  ret(0) = tl;
  ret(1) = er;
  ret(2) = sm;
  return(ret);
}



// [[Rcpp::export]]
int select_genes( arma::mat& data, int row, int col, int g, int k, int r, double b1, int b2)
{
  //  ("How many random starts for the fitting of t components to individual genes?\n"); r
  // ("How many k-means starts for the fitting of t components to individual genes?\n"); k
  //printf ("Enter threshold for likelihood ratio statistic:\n");	/* -2 log \lambda cutoff */ b1
  // printf ("Enter threshold for minimum cluster size:\n");	/* s_{min} cutoff */ b2

  int x;			/* counter */
  int i;			/* another counter */
  int rc, r1, k1, clus;
  int genecount = 0;

   arma::mat  sg =  arma::zeros<arma::mat>(MAX_GENES,col);
   arma::mat  f4 =  arma::zeros<arma::mat>(MAX_GENES,col);
   arma::mat  f3 =  arma::zeros<arma::mat>(row,3);

  //struct tab gtab[MAX_SELECTED_GENES];	/* for our cut-down genes and -2 log \lambdas */
  std::vector<tab> gtab;
  g = 1; //why

  for(int rc=0;rc<row;rc++)
  {
     arma::vec t =   arma::zeros(col);
    double last;
    int sm, er;
    double tl;

    t = data.row(rc);

    // int row;			/* 2nd argument */ argv
    // int col;			/* 3rd argument */
    // int g;			/* 4th argument */

     arma::vec em_res = call_emmix_sel(t, row, col, g);	/* get -2 log \lambda, smaller, error */ //call emmix
    tl = em_res(0);
    er = em_res(1);
    sm = em_res(2);

    if (tl > b1)		/* if -2 log \lambda exceeds b_1 cutoff */
    {
      if (sm >= b2)		/* if smaller groups exceeds b_2 cutoff */
      {
         arma::vec temp =  arma::zeros(3);
        temp(0) = rc;
        temp(1) = tl;
        temp(2) = sm;
        f3.row(rc) = temp;

        gtab[genecount].nu = genecount; //results
        gtab[genecount].tl = tl;
        gtab[genecount].sm = sm;

        sg.row(genecount)= t; //copies s to sg[genecount]
        genecount++;

      } else {

        /* g is the number of groups */
         arma::vec em_res = call_emmix_sel(t, row, col, g+1); //same but g+1
        tl = em_res(0);
        er = em_res(1);
        sm = em_res(2);

        if (tl > b1)
        {
          arma::vec temp =  arma::zeros(3);
          temp(0) = rc;
          temp(1) = tl;
          temp(2) = sm;
          f3.row(rc) = temp;

          gtab[genecount].nu = genecount;
          gtab[genecount].tl = tl;
          gtab[genecount].sm = sm;
          sg.row(genecount)= t;
          genecount++;
        }
      }

    }

  }


  std::sort(gtab.begin(), gtab.end(), tabcompare);


  // save data
  arma::vec tmp = arma::zeros(3);

  for(int j=0;j<genecount;j++){
    f4.row(j)=sg.row(gtab[j].nu);
    tmp(0)=gtab[j].nu+1;
    tmp(1)=gtab[j].tl;
    tmp(2)=gtab[j].sm;
    f3.row(j)=tmp;
  }

  //need to write or save f3 and f4

  return 0;
}



