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
#' @param alpha the significance level for accuracy measures. Default value is 0.05
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
#' @include utils.R
#'
#' @import HDInterval
#'
#' @export
postNtcondN0 <- function(N0, nMNOmat,  distNames, variation, n = 1e3, alpha = 0.05) {

  Ntmat <- rNtcondN0(n, N0, nMNOmat, distNames, variation)
  postMean <- apply(Ntmat, 2, function(N){round(mean(N))})
  postSD <- apply(Ntmat, 2, function(N){round(sd(N))})
  postCV <- apply(Ntmat, 2, function(N){round(sd(N)/mean(N), 2)})


  postMedian <- apply(Ntmat, 2, function(N){round(median(N))})
  postMedian_CILB <- apply(Ntmat, 2, function(N){equalTailedInt(N, alpha)['lower']})
  postMedian_CIUB <- apply(Ntmat, 2, function(N){equalTailedInt(N, alpha)['upper']})
  postMedianQuantileCV <- apply(Ntmat, 2, function(N){round( IQR(N)/median(N) *100 ,2)})

  postMode <- apply(Ntmat, 2, function(N){Mode(N)})

  postMode_CILB <- apply(Ntmat, 2, function(N){hdi(N, 1-alpha)['lower']})
  postMode_CIUB <- apply(Ntmat, 2, function(N){hdi(N, 1-alpha)['upper']})
  postModeQuantileCV <- apply(Ntmat, 2, function(N){round(IQR(N)/Mode(N) *100,2)})

  output <- list(postMean = postMean, postSD = postSD, postCV = postCV, postMedian = postMedian, postMedian_CILB = postMedian_CILB, postMedian_CIUB = postMedian_CIUB, postMedianQuantileCV = postMedianQuantileCV, postMode = postMode, postMode_CILB = postMode_CILB, postMode_CIUB = postMode_CIUB, postModeQuantileCV = postModeQuantileCV)
  output <- Reduce(cbind, output)
  colnames(output) <- c('postMean', 'postSD', 'postCV', 'postMedian', 'postMedian_CILB', 'postMedian_CIUB', 'postMedianQuantileCV', 'postMode', 'postMode_CILB', 'postMode_CIUB', 'postModeQuantileCV')
  return(output)

  return(postMean)
}

