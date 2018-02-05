#' @title kummmer
#' @description a function
#' @author David Salgado, Bogdan Oancea
#' @export
#'

kummer <- function(x, a, b, relTol = 1e-6){

  output <- Kummer(x, a, b, relTol)
  output[is.infinite(output)] <- .Machine$double.xmax
  return(output)
}
