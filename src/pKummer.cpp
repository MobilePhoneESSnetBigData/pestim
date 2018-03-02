// [[Rcpp::depends(RcppParallel)]]
#include "RcppParallel.h"
#include "Rcpp.h"
#include <cmath>
using namespace RcppParallel;



struct wKummer : public Worker {
  
  //input
  const RVector<double> z;
  const RVector<double> a;
  const RVector<double> b;
  const double relTol;
  
  // output
  RVector<double> suma;
  
  // initialize from Rcpp input and output
  wKummer(Rcpp::NumericVector z, Rcpp::NumericVector a, Rcpp::NumericVector b, double relTol, Rcpp::NumericVector suma)
    : z(z), a(a), b(b), relTol(relTol), suma(suma) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      suma[i] = 1;
      double sumando = 1;
      //double dif = a[i] - b[i];
      //double dif2 = 1- a[i];
      
      if (z[i]<80) {
        for (long j=0;std::fabs(sumando/suma[i])>relTol;++j){
          sumando*=(z[i] / (j+1)) * ((a[i] + j) / (b[i] + j));
          suma[i]+=sumando;
        }
      }
      else {
        for (long j=0;fabs(sumando/suma[i])>relTol;++j){
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
  
  wKummer wk(z, a, b, relTol, suma);
  
  // call it with parallelFor
  parallelFor(0, z.length(), wk);
  
  return(suma);
}
