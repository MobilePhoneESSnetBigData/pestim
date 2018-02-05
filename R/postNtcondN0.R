#' @title Time evolution of the posterior mean, median, and mode for the number of individuals conditional on their estimated initial values
#' @description Compute the time evolution of the posterior mean, median, and mode for the number of individuals
#' conditional on their estimated initial values, according to the hierarchical model.
#' @author David Salgado
#'
#' @param N0 initial number of population in each cell
#'
#' @param nMNOmat
#'
#' @param distNames
#'
#' @param variation
#'
#' @param n number of points to generate in the posterior distribution for the computation. Default
#' value is 1e3
#'
#' @return
#'
#' @include rNtcondN0.R
#'
#' @export
#'

postNtcondN0 <- function(N0, nMNOmat,  distNames, variation, n = 1e3){

  Ntmat <- rNtcondN0(n, N0, nMNOmat, distNames, variation)
  postMean <- apply(Ntmat, 2, function(N){round(mean(N))})
  postMedian <- apply(Ntmat, 2, function(N){round(median(N))})
  postMode <- apply(Ntmat, 2, function(N){N[which.max(names(table(N)))]})
  output <- list(postMean = postMean, postMedian = postMedian, postMode = postMode)
  output <- Reduce(cbind, output)
  colnames(output) <- c('postMean', 'postMedian', 'postMode')
  return(output)

  return(postMean)
}
