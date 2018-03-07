#include <Rcpp.h>
#include <cmath>
#include "compute.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Kummer(NumericVector z, NumericVector a, NumericVector b, double relTol) {
  //double sumando;
  int n=z.size();
  //long j;
  NumericVector suma(n);
  compute(z,a,b,suma, relTol, 0, n);
  // for (i=0;i<n;++i) {
  //   if (z[i]<80) {
  //     for (j=0,sumando=1,suma[i]=1;fabs(sumando/suma[i])>relTol;++j){
  //       sumando*=(z[i] / (j+1)) * ((a[i] + j) / (b[i] + j));
  //       suma[i]+=sumando;
  //     }
  //   }
  //   else {
  //     for (j=0,sumando=1,suma[i]=1;fabs(sumando/suma[i])>relTol;++j){
  //       sumando*=((1-a[i]+j) / (j+1)) * ((b[i] - a[i] + j) / (z[i]));
  //       suma[i]+=sumando;
  //     }
  //     suma[i]*=exp(z[i]+(a[i]-b[i])*log(z[i])+lgamma(b[i])-lgamma(a[i]));
  //   }
  // }

  return(suma);
}
