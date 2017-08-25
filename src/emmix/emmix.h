#ifndef EMMIX_H
#define EMMIX_H

#ifdef __cplusplus
extern "C" {
#endif

///expose only fortran fuctions we want

  // SUBROUTINE MODEL
  // Obtain parameters defining crustal motion model
  void model_();

  // SUBROUTINE NEWCOR
  void main_(NIND,NATT,NCOV,X,NG0,NG1,XMU,XVAR,T,MAXITS, TOLS,AIT,OPTION,START,NRANDS,NKMEANS,FACT,RDSB,NATTQ);

  //void iymdmj_(const int *iyr, const int *imon, const int *iday, int *mjd);

  //how to figure out what var types we have, also everything is a pointer
  SUBROUTINE FIT(NIND,NATT,NCOV,X,NG,XMU,XVAR,T,XMU2,XVAR2,T2, TOL,MAXIT,AIT,XLIKE,TLIKE,TMP,FACT,B,D,NATTQ,IER)

#ifdef __cplusplus
}
#endif

#endif // EMMIX_H
