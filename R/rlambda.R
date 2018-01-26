#' @title rlambda
#' @description
#' @author David Salgado
#' @export
#'
rlambda <- function(n, nMNO, nReg, fu, fv, flambda, relTol = 1e-2, nSim = 1e3, nStrata = c(1, 1e2), verbose = FALSE){

  nCells <- length(nMNO)
  if (nCells != 1) stop('Only one cell at a time.')
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')

  ######  Computing betaOpt
  lambdaOpt <- modeLambda(nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose)

  ###### Computing rejection rate  #####################
  if (verbose) cat('Computing rejection rate...')
  f <- function(x){dlambda(x, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose)$probLambda}
  location <- lambdaOpt
  scale <- sqrt(flambda$shape * flambda$scale^2)
  F0 <- pcauchy(0, location = location, scale = scale)
  g <- function(x){
    #dgamma(x, shape = nMNO + 1, scale = lambdaOpt / nMNO)
    dcauchy(x, location = location, scale = scale) / (1 - F0)
    #dg(x, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata)
  }
  fun <- function(x){g(x) / f(x)}

  optimC <- 1.05 / optimise(fun, interval = c(max(location - scale, 0), location + scale))$objective
  if (verbose) cat(' ok.\n')

  if (verbose) cat('Generating and accepting/rejecting values...\n')
  u <- runif(2 * n)
  x <- qcauchy(F0 + u * ( 1 - F0), location = location, scale = scale)
  #x <- rgamma(10 * n, shape = nMNO + 1, scale = lambdaOpt / nMNO)
  v <- runif(n)
  if (verbose) cat('   of target distribution...')
  fx <- f(x)
  if (verbose) cat(' ok.\n')
  if (verbose) cat('   of candidate distribution...')
  gx <- dcauchy(x, location = location, scale = scale)
  if (verbose) cat(' ok.\n')
  output <- x[v <= fx / (optimC * gx)]
  if (verbose) cat(paste0(length(output), ' points selected.\n'))
  while (length(output) < n) {
    u <- runif(2 * n)
    x <- qcauchy(F0 + u * ( 1 - F0), location = location, scale = scale)
    #x <- rgamma(10 * n, shape = nMNO + 1, scale = lambdaOpt / nMNO)
    v <- runif(n)
    fx <- f(x)
    gx <- g(x)
    aux <- x[v <= fx / (optimC * gx)]
    if (verbose) cat(paste0(length(output), ' points selected.\n'))
    output <- c(output, aux)
  }
  output <- output[1:n]
  if (verbose) cat(' ok.\n')
  return(output)
}
