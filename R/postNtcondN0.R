#' @title postNtcondN0
#' @description postNtcondN0
#' @author David Salgado
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
