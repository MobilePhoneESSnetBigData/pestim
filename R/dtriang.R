#' @title dtriang
#' @description
#' @author David Salgado
#' @export
#'
dtriang <- function(x, xMin, xMax, xMode){

  n <- length(x)
  if (length(xMin) == 1) xMin <- rep(xMin, n)
  if (length(xMax) == 1) xMax <- rep(xMax, n)
  if (length(xMode) == 1) xMode <- rep(xMode, n)

  output <- x
  output[x <= xMin] <- 0
  output[x >= xMax] <- 0
  range1 <- (x > xMin & x <= xMode)
  output[range1] <- (2 * (x[range1] - xMin[range1])) / ((xMax[range1] - xMin[range1]) * (xMode[range1] - xMin[range1]))
  range2 <- (x >= xMode & x < xMax)
  output[range2] <- (2 * (xMax[range2] - x[range2])) / ((xMax[range2] - xMin[range2]) * (xMax[range2] - xMode[range2]))
  return(output)
}
