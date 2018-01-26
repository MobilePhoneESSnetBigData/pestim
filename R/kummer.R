#' @title kummmer
#' @description
#' @author David Salgado
#' @export
#'
Rcpp::cppFunction('
                  double atomKummer(double x, double a, double b, double relTol) {
                  double sumando, suma;
                  int n;

                  for(n=1, sumando=1, suma=1; sumando / suma > relTol; ++n, sumando*= x * (a+n-1) / ((b+n-1)*n), suma+=sumando);
                  return(suma);
                  }
                  ')

kummer <- function(x, a, b, relTol = 1e-6){

  output <- sapply(seq(along = x), function(i){atomKummer(x[i], a[i], b[i], relTol)})
  output[is.infinite(output)] <- .Machine$double.xmax
  return(output)

}
