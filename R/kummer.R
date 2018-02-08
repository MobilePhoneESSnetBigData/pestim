#' @title kummmer
#' @description
#' @author David Salgado, Bogdan Oancea
#' @useDynLib pestim
#' @export
#'
kummer <- function(x, a, b, relTol = 1e-6){

  output <- atomKummer(x, a, b, relTol)
  output[is.infinite(output)] <- .Machine$double.xmax
  return(output)

}
