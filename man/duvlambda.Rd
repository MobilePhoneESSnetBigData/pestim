% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duvlambda.R
\name{duvlambda}
\alias{duvlambda}
\title{Posterior density function of the parameters u, v, and lambda.}
\usage{
duvlambda(lambda, u, v, nMNO, fu, fv, flambda, relTol = 1e-06,
  nThreads = RcppParallel::defaultNumThreads())
}
\arguments{
\item{lambda}{values of the parameter \eqn{\lambda}}

\item{u}{values of the parameter \eqn{u}}

\item{v}{values of the parameter \eqn{v}}

\item{nMNO}{non-negative integer vectors with the number of individuals detected according to the
network operator}

\item{fu}{named list with the prior marginal distributions of the hyperparameter \eqn{u}}

\item{fv}{named list with the prior marginal distributions of the hyperparameter \eqn{v}}

\item{flambda}{named list with the prior marginal distribution of the parameter \eqn{\lambda}}

\item{relTol}{relative tolerance in the computation of the \code{\link{kummer}} function. Default
value is \code{1e-6}}

\item{nThreads}{number (default the number of all cores, including logical cores) to use for 
computation}
}
\value{
\code{duvlambda} returns the probability mass function of the number of indviduals 
detected by the mobile network operator according to the hierarchical model. It depends on priors
for the parameters \deqn{u}, \deqn{v}, \deqn{\lambda}.
}
\description{
Compute the unnormalized posterior density function of the parameter \eqn{\lambda}
in the hierarchical model.
}
\details{
The prior distributions are specified as named lists where the first component of each 
list must be the name of distribution ('unif', 'triang', 'gamma') and the rest components must be
named according to the name of the parameters of the density/probability function of the 
corresponding distribution according to:

  \itemize{

    \item unif: \code{xMin}, \code{xMax} for the minimum, maximum of the sampled interval.
    \item triang: \code{xMin}, \code{xMax}, \code{xMode} for minimum, maximum and mode (see
    \code{\link{qtriang}}).
    \item gamma: \code{scale} and \code{shape} with the same meaning as in \code{\link{rgamma}}.
  }
}
\examples{
f <- function(x){duvlambda(x, 0.346, 97, 19, 
           fu = list('unif', xMin = 0.3, xMax = 0.4),
           fv = list('unif', xMin = 90, xMax = 105),
           flambda = list('gamma', shape = 11, scale = 97 / 10))$prob}
curve(f, 0, 150)

f <- function(x){duvlambda(70, x, 97, 19, 
           fu = list('unif', xMin = 0.3, xMax = 0.4),
           fv = list('unif', xMin = 90, xMax = 105),
           flambda = list('gamma', shape = 11, scale = 97 / 10))$prob}
curve(f, 0, 1)

f <- function(x){duvlambda(70, 0.35, x, 19, 
           fu = list('unif', xMin = 0.3, xMax = 0.4),
           fv = list('unif', xMin = 90, xMax = 105),
           flambda = list('gamma', shape = 11, scale = 97 / 10))$prob}
curve(f, 80, 115)

}
\seealso{
\code{\link{dlambda}}
}
