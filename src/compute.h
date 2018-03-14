#ifndef __COMPUTE_H__
#define __COMPUTE_H__
#include <cmath>

template<typename T>
  void compute(const T& z, const T& a, const T& b, T& suma, double relTol, std::size_t begin, std::size_t end) {
    double sumando;
    std::size_t i;
    long j;
    for (i=begin;i<end;++i) {
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
#endif
