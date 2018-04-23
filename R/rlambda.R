#' @title Generation of random deviates of the posterior distribution of parameter lambda.
#'
#' @description Generate random points according to the posterior probability distribution of the
#' parameter lambda in the hierarchical model.
#'
#' @param n number of values to generate
#'
#' @param nMNO, nReg non-negative integer vectors with the number of individuals detected in each
#' cell according to the network operator and the register
#'
#' @param fu, fv named lists with the prior marginal distributions of the two-dimensional points
#' for the Monte Carlo integration
#'
#' @param flambda named list with the prior distribution of the lambda parameter
#'
#' @param relTol relative tolerance in the computation of the \code{\link{kummer}} function. Default
#' value is \code{1e-6}
#'
#' @param nSim number of two-dimensional points to generate to compute the integral. Default value
#' is \code{1e4}
#'
#' @param nStrata integer vector of length 2 with the number of strata in each dimension. Default
#' value is \code{c(1, 1e2)}
#'
#' @param verbose logical (default \code{FALSE}) to report progress of the computation
#'
#' @return \code{rlambda} generates \code{n} points according to the posterior distribution of
#' the parameter lambda. The function returns a vector with these points.
#'
#' @details The points are generated according to the accept-reject method using as candidate
#' distribution a Cauchy distribution whose parameters are taken from the prior distributions and
#' the mode of the posterior distribution of the lambda parameter.
#'
#' The prior distributions are specified as named lists where the first component of each list must
#' be the name of distribution ('unif', 'triang', 'degen', 'gamma') and the rest components must be
#' named according to the name of the parameters of the random generator of the corresponding
#' distribution according to:
#'
#'   \itemize{
#'
#'     \item unif: \code{xMin}, \code{xMax} for the minimum, maximum of the sampled interval.
#'     \item degen: \code{x0} for the degenerate value of the random variable.
#'     \item triang: \code{xMin}, \code{xMax}, \code{xMode} for minimum, maximum and mode (see
#'     \code{\link{qtriang}}).
#'     \item gamma: \code{scale} and \code{shape} with the same meaning as in \code{\link{rgamma}}.
#'   }

#'
#' @seealso \code{\link{dlambda}}, \code{\link{rg}} for related functions.
#'
#' @examples
#' # It takes a couple of minutes
#' hist(rN0(500, nMNO = 20, nReg = 115, fu = list('unif', xMin = 0.3, xMax = 0.5),
#'         fv = list('unif', xMin = 100, xMax = 120),
#'         flambda = list('gamma', shape = 11, scale = 12))$N0,
#'         breaks = seq(1, 200, by = 1), main ='', xlab = 'number of individuals')
#'
#' @include modeLambda.R
#'
#' @export
rlambda <- function(n, nMNO, nReg, fu, fv, flambda, relTol = 1e-6, nSim = 1e4, nStrata = c(1, 1e2), verbose = FALSE, nThreads = RcppParallel::defaultNumThreads()){

  nCells <- length(nMNO)
  if (length(nReg) != nCells) stop('nReg and nMNO must have the same length.')
  #if (nCells != 1) stop('Only one cell at a time.')

  if (nCells == 1){

    ######  Computing lambdaOpt
    lambdaOpt <- modeLambda(nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)

    ###### Computing rejection rate  #####################
    if (verbose) cat('Computing rejection rate...')
    f <- function(x){dlambda(x, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata, verbose, nThreads)$probLambda}

    if (flambda[[1]] == 'gamma') {

      varfLambda <- (flambda[['shape']] - 1) * flambda[['scale']]**2
      alpha_candidate <- lambdaOpt**2 / varfLambda

    }

    shape_candidate <- alpha_candidate + 1
    scale_candidate <- lambdaOpt / alpha_candidate

    location <- lambdaOpt
    scale <- sqrt(flambda$shape * flambda$scale^2)
    g <- function(x){
      dgamma(x, shape = shape_candidate, scale = scale_candidate)
      #dcauchy(x, location = location, scale = scale) / (1 - F0)
      #dg(x, nMNO, nReg, fu, fv, flambda, relTol, nSim, nStrata)

    }

    fun <- function(x){g(x) / f(x)}
    optimC <- f(lambdaOpt) / g(lambdaOpt)

    #optimC <- 1 / optimise(fun, interval = c(max(location - scale, 0), location + scale))$objective
    if (verbose) cat(' ok.\n')

    nEffect <- ceiling(1.5 * n / (1 - optimC)) # Factor 1.5 to overproduce candidate points

    if (verbose) cat('Generating and accepting/rejecting values...\n')
    #u <- runif(2 * n)
    #x <- qcauchy(F0 + u * ( 1 - F0), location = location, scale = scale)
    x <- rgamma(nEffect, shape = shape_candidate, scale = scale_candidate)
    v <- runif(nEffect)
    if (verbose) cat('   of target distribution...\n')
    fx <- f(x)
    if (verbose) cat('   ok.\n')
    if (verbose) cat('   of candidate distribution...\n')
    #gx <- dcauchy(x, location = location, scale = scale)
    gx <- dgamma(x, shape = shape_candidate, scale = scale_candidate)
    if (verbose) cat('   ok.\n')
    output <- x[v <= fx / (optimC * gx)]
    if (verbose) cat(paste0(length(output), ' points selected.\n'))
    while (length(output) < n) {

      nEffect <- ceiling(1.5 * (n - length(output)) / (1 - optimC)) # Factor 1.5 to overproduce candidate points

      #u <- runif(2 * n)
      #x <- qcauchy(F0 + u * ( 1 - F0), location = location, scale = scale)
      x <- rgamma(nEffect, shape = shape_candidate, scale = scale_candidate)
      v <- runif(nEffect)
      fx <- f(x)
      gx <- g(x)
      aux <- x[v <= fx / (optimC * gx)]
      if (verbose) cat(paste0(length(output), ' points selected.\n'))
      output <- c(output, aux)
    }
    output <- output[1:n]
    if (verbose) cat(paste0(length(output), ' points selected.\n'))
    if (verbose) cat(' ok.\n')
    return(output)

  } else {

    output <- lapply(seq(along = nMNO), function(i){

      rlambda(n, nMNO[i], nReg[i], fu[[i]], fv[[i]], flambda[[i]], relTol, nSim, nStrata, verbose, nThreads)

    })

    output <- Reduce(cbind, output)
    dimnames(output) <- NULL
    return(output)
  }
}
