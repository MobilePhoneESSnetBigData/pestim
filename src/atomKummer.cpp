#include<Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
                  double atomKummer(double x, double a, double b, double relTol) {
                  double sumando, suma;
                  int n;

                  for(n=1, sumando=1, suma=1; sumando / suma > relTol; ++n, sumando*= x * (a+n-1) / ((b+n-1)*n), suma+=sumando);
                  return(suma);
                  }
