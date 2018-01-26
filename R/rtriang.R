#' @title rtriang
#' @description
#' @author David Salgado
#' @export
#'
rtriang <- function(n, xMin, xMax, xMode){

  u <- runif(n)
  mc <- match.call()
  mc[[1L]] <- qtriang
  mc[['n']] <- NULL
  mc[['u']] <- u
  output <- eval(mc)
  return(output)
}
