#' @title kummmer
#' @description
#' @author David Salgado, Bogdan Oancea
#' @export
#'
kummer2 <- function(x, a, b, relTol = 1e-6){

  output <- sapply(seq(along = x), function(i){atomKummer(x[i], a[i], b[i], relTol)})
  output[is.infinite(output)] <- .Machine$double.xmax
  return(output)

}

kummer <- function(x, a, b, relTol = 1e-6){
  
  output <- atomKummer(x, a, b, relTol)
  output[is.infinite(output)] <- .Machine$double.xmax
  return(output)
  
}