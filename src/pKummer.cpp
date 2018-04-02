#include "Rcpp.h"
#include "RcppParallel.h"
#include "compute.h"

using namespace RcppParallel;


// [[Rcpp::depends(RcppParallel)]]
struct wKummer : public Worker {

  //input
  const RVector<double> z;
  const RVector<double> a;
  const RVector<double> b;
  const double relTol;

  // output
  RVector<double> suma;
  RVector<int> err;

  // initialize from Rcpp input and output
  wKummer(const Rcpp::NumericVector& z1, const Rcpp::NumericVector& a1, const Rcpp::NumericVector& b1, const double relTol1, Rcpp::NumericVector& suma1, Rcpp::IntegerVector& error)
    : z(z1), a(a1), b(b1), relTol(relTol1), suma(suma1), err{error} {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    compute<RVector<double>, RVector<int>>(z,a,b,suma, relTol, begin, end, err);
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector pKummer(const Rcpp::NumericVector& z, const Rcpp::NumericVector& a, const Rcpp::NumericVector& b, double relTol, Rcpp::IntegerVector& err) {
  int n=z.size();
  Rcpp::NumericVector suma(n);
  int grain_size = 1;

  // an automatic procedure to compute the optimum grain size  - to be developed
  // for now we use some heuristic numbers
  if (n >= 1e4)
    grain_size = 2000;
  else
    grain_size = 500;


  wKummer wk(z, a, b, relTol, suma, err);

  // call it with parallelFor
  parallelFor(0, z.length(), wk, grain_size);

  return(suma);
}
