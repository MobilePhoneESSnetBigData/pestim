#ifndef __COMPUTE_H__
#define __COMPUTE_H__
#include <cmath>
#include <climits>

#define SUCCESS 0
#define NO_CONVERGENCE 1
#define OVERFLOW_ERROR 2

template<typename T, typename S>
  void compute(const T& z, const T& a, const T& b, T& suma, double relTol, std::size_t begin, std::size_t end, S& error) {
    double a1;
    std::size_t i;
    int j;
    for ( i = begin; i < end; ++i) {
      a1 = 1;
      suma[i]=1;
      for(j = 0; j < 1000; ++j) {
        a1 *= (a[i]+j)/(b[i]+j) * z[i]/(j+1);
        suma[i] += a1;
        if (fabs(a1)/fabs(suma[i]) < relTol )
          break;
      }
      if(j == 1000) {
        error[i] = NO_CONVERGENCE;
      }
      if(suma[i] == INFINITY || suma[i] == NAN)
        error[i] = OVERFLOW_ERROR;
    }
  }
#endif
