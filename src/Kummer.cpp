#include <Rcpp.h>
#include <cmath>
#include "compute.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Kummer(const NumericVector& z, const NumericVector& a, const NumericVector& b, const double relTol) {
  int n=z.size();
  NumericVector suma(n);
  compute<NumericVector>(z,a,b,suma, relTol, 0, n);
  return(suma);
}
