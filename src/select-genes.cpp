//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//'@importFrom Rcpp sourceCpp
//'@useDynLib EMMIXgene

#include <vector>
#include <sting>


#define LINE_LENGTH 5000	/* max length of each line in microarray data */
#define MAX_GENES 5000		/* max number of genes in microarray data */
#define MAX_TISSUES 200     /* max number of tissues in microarray data */
#define MAX_GROUPS 100		/* max number of groups to split into */
#define MAX_ENTRY_LENGTH 20	/* max length of any single entry in the array */


struct tab
{
  int nu;
  int count;
  double tl;
  std::vector<double> s;
};



// [[Rcpp::export]]
int select_genes(std::string filename, int row, int col, int g, int k, int r, double b1, int b2)
{
  //  ("How many random starts for the fitting of t components to individual genes?\n"); r
  // ("How many k-means starts for the fitting of t components to individual genes?\n"); k
  //printf ("Enter threshold for likelihood ratio statistic:\n");	/* -2 log \lambda cutoff */ b1
  // printf ("Enter threshold for minimum cluster size:\n");	/* s_{min} cutoff */ b2

  int x;			/* counter */
  int i;			/* another counter */
  int rc, b1, b2, r1, k1, clus;
  int genecount = 0;
  arma::vec s =  arma::zeros(LINE_LENGTH);
  arma::vec sg = arma::zeros(MAX_GENES);

  struct tab gtab[MAX_SELECTED_GENES];	/* for our cut-down genes and -2 log \lambdas */
  g = 1; //why

  rc = 0;
  while (rc < row)
  {
    arma::vec t =  arma::zeros(LINE_LENGTH);
    double last;
    int sm, er;
    double tl;

    /* grab line i, change space to CR, remove blank lines */
    /* sed -n $i'p' $1 | tr ' ' '\n' | sed '/^$/d' > tmp.$$ */
    /* write output to tmp.$$ */

    fgets (s, LINE_LENGTH, f); //grab row - how to ensue is the right one
    rc++;

    call_emmix (argv, t, &tl, &sm, &er, col, g);	/* get -2 log \lambda, smaller, error */ //call emmix

    if (tl > b1)		/* if -2 log \lambda exceeds b_1 cutoff */
    {
      if (sm >= b2)		/* if smaller groups exceeds b_2 cutoff */
      {
        fprintf (f3, "%d %f %d\n", rc, tl, sm); //saves result

        gtab[genecount].nu = genecount; //results
        gtab[genecount].tl = tl;
        gtab[genecount].sm = sm;

        strcpy (sg[genecount], s); //copies s to sg[genecount]
        genecount++;

      } else {

        /* g is the number of groups */
        call_emmix (argv, t, &tl, &sm, &er, col, g + 1); //same but g+1
        if (tl > b1)
        {
          fprintf (f3, "%d %f %d\n", rc, tl, sm);

          gtab[genecount].nu = genecount;
          gtab[genecount].tl = tl;
          gtab[genecount].sm = sm;
          strcpy (sg[genecount], s);
          genecount++;
        }
      }

    }

  }

  qsort ((void *) gtab, genecount, sizeof (struct tab), tabcompare); //sort gtab on .tl


  // save data
  rc = 0;
  while (rc < genecount){
    fputs (sg[gtab[rc++].nu], f4);
  }

  rc = 0;
  while (rc < genecount)
  {
    fprintf (f3,"%d %f %d\n",gtab[rc].nu+1, gtab[rc].tl, gtab[rc].sm);
    rc++;
  }


  return 0;
}


