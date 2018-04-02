#include <Rcpp.h>
#include <cmath>
#include "compute.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Kummer(const NumericVector& z, const NumericVector& a, const NumericVector& b, const double relTol, NumericVector& err) {
  int n=z.size();
  NumericVector suma(n);
  compute<NumericVector, NumericVector>(z,a,b,suma, relTol, 0, n, err);
  return(suma);
}
