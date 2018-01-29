#include<Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
NumericVector atomKummer(NumericVector x, NumericVector a, NumericVector b, double relTol) {
  int n;
  int l = x.size();
  NumericVector sumando(l), suma(l);

  long j;
  for(j=0;j<l;j++)
    for(n=1, sumando[j]=1, suma[j]=1; sumando[j] / suma[j] > relTol; ++n, sumando[j] *= x[j] * (a[j]+n-1) / ((b[j]+n-1)*n), suma[j] +=sumando[j]);

  return(suma);
}
