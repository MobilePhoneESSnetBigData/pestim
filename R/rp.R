#' @title kummmer
#' @description postestimates of
#' @author David Salgado
#' @export
#'

rp <- function(n, flist){

  alphaParam <- genAlpha(n, flist)
  output <- t(apply(alphaParam, 1, function(alphaPar){MCMCpack::rdirichlet(n = 1, alpha = alphaPar)}))
  return(output)

}
