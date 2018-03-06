#include "Rcpp.h"
#include "RcppParallel.h"
#include <cmath>
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

  // initialize from Rcpp input and output
  wKummer(Rcpp::NumericVector z1, Rcpp::NumericVector a1, Rcpp::NumericVector b1, double relTol1, Rcpp::NumericVector suma1)
    : z(z1), a(a1), b(b1), relTol(relTol1), suma(suma1) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    int j;
    double sumando;
    for (std::size_t i = begin; i < end; i++) {
      if (z[i]<80) {
        for (j=0,sumando=1,suma[i]=1;fabs(sumando/suma[i])>relTol;++j){
          sumando*=(z[i] / (j+1)) * ((a[i] + j) / (b[i] + j));
          suma[i]+=sumando;
        }
      }
      else {
        for (j=0,sumando=1,suma[i]=1;fabs(sumando/suma[i])>relTol;++j){
          sumando*=((1-a[i]+j) / (j+1)) * ((b[i] - a[i] + j) / (z[i]));
          suma[i]+=sumando;
        }
        suma[i]*=exp(z[i]+(a[i]-b[i])*log(z[i])+lgamma(b[i])-lgamma(a[i]));
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector pKummer(Rcpp::NumericVector z, Rcpp::NumericVector a, Rcpp::NumericVector b, double relTol) {
  int n=z.size();
  Rcpp::NumericVector suma(n);
  int grain_size = 1;

  // an automatic procedure to compute the optimum grain size  - to be developed
  // for now we use some heuristic numbers
  if (n >= 1e4)
    grain_size = 1000;
  else
    grain_size = 100;


  wKummer wk(z, a, b, relTol, suma);

  // call it with parallelFor
  parallelFor(0, z.length(), wk, grain_size);

  return(suma);
}
