#' @title Utility function used throughout the package.
#'
#' @description Containts functions to compute the Mode of a data set and the equal Tailed interval
#'
#' @export
Mode = function(v){
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

equalTailedInt <- function(x, alpha){
  output <- quantile(x, c((1 - alpha) / 2, (1 + alpha) / 2))
  names(output) <- c('lower', 'upper')
  return(output)
}
