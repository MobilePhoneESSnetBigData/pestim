#include<Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
NumericVector atomKummer(NumericVector x, NumericVector a, NumericVector b, double relTol) {
  long n;
  long dim = x.size();
  NumericVector sumando(dim,1.), suma(dim,1.);

  long j;
  for(j=0;j<dim;j++) {
    a[j]--;
    b[j]--;
    for(n=1; sumando[j] / suma[j] > relTol; ++n, sumando[j] *= x[j] * (a[j]+n) / ((b[j]+n)*n), suma[j] +=sumando[j]);
  }
  return(suma);
}
