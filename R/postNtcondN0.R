#' @title Posterior mean, median, and mode for the number of individuals at an arbitrary time
#' conditioned upon the initial population.
#'
#' @description Compute the posterior mean, median, and mode for the number of individuals
#' generating posterior distribution according to the hierarchical model conditioned upon the
#' initial population of each cell, which must be provided
#'
#' @param N0 initial population in each cell
#'
#' @param nMNOmat transition matrix with the number of individuals displaced from cell to cell
#' detected by the Mobile Network Operator
#'
#' @param distNames character vector with the names of the prior distributions for each cell
#'
#' @param variation list of lists whose components are parameters providing a measure of variation
#' of each prior distribution
#'
#' @param n number of points to generate in the posterior distribution for the computation. Default
#' value is 1e3
#'
#' @return Return a matrix with three columns (mean, median, and mode estimates) and one row per
#' cell
#'
#' @examples
#' ## First, the inputs:
#'
#' # The initial population
#' N0 <- c(93, 123, 130)
#'
#' #The transition matrix of individuals detected by the MNO
#' nMNOmat <- rbind(c(10, 3, 4), c(5, 21, 3), c(3, 9, 18))
#'
#' # Names and parameters of priors for the transition probabilities
#' distNames <- rep('unif', 3)
#' variation <- rep(list(list(cv = 0.20)), 3)
#'
#' # It takes a couple of minutes.
#' postNtcondN0(N0, nMNOmat, distNames, variation)
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
