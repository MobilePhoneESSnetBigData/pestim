#' @title qtriang
#' @description
#' @author David Salgado
#' @export
#'
qtriang <- function(u, xMin, xMax, xMode){

  n <- length(u)
  if (length(xMin) == 1) xMin <- rep(xMin, n)
  if (length(xMax) == 1) xMax <- rep(xMax, n)
  if (length(xMode) == 1) xMode <- rep(xMode, n)

  output <- u
  range1 <- (u < (xMode - xMin) / (xMax - xMin))
  output[range1] <- xMin[range1] + sqrt(u[range1] * (xMax[range1] - xMin[range1]) * (xMode[range1] - xMin[range1]))
  range2 <- ( u > (xMode - xMin) / (xMax - xMin))
  output[range2] <- xMax[range2] - sqrt((1 - u[range2]) * (xMax[range2] - xMin[range2]) * (xMax[range2] - xMode[range2]))
  return(output)
}
