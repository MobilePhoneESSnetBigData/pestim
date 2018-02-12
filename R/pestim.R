#' pestim: a hierarchical model to estimate population counts with aggregated mobile phone data.
#'
#' This package provides an implementation for a hierarchical model to combine both aggregated
#' mobile phone data and external official (administrative or survey) data to produce estimates
#' of population counts in each cell of a division of a territory.
#'
#' @section Context:
#'
#' This package has been developed in the context of a European research project within the European
#' Statistical System called
#' \href{https://webgate.ec.europa.eu/fpfis/mwikis/essnetbigdata/index.php/Main_Page}{ESSnet on Big
#' Data}. More specifically this work corresponds to the work package on mobile phone data by which
#' we assess the use of this data source in the production of official statistics. The goals of the
#' project is many-fold. Firstly, the issue of accessing these data for the production of official
#' statistics initially for research and then for standard production has been investigated.
#' Secondly, in a hands-on bottom-up approach, we make some initial methodological proposals to
#' produce concrete statistical output using those data sets compiled in the preceding phase.
#' Thirdly, in parallel, IT tools, architecture and software development are assessed especially in
#' contrast to traditional computer frameworks. Finally, quality is appraised especially in the
#' context of the European Statistics Code of Practice and ESS Quality Assurance Framework. This
#' package provides a first-step implementation of software routines to present a proof of concept
#' about a methodological proposal (see below) to make inferences about a target population from a
#' mobile phone dataset.
#'
#' @section The hierarchical model in a nutshell:
#' The methodological proposal giving rise to this package focuses on the inference exercise
#' connecting aggregated mobile phone data with a target population under analysis. In concrete, the
#' goal is to provide estimates of population counts in each cell in which we have divided the
#' territory for which the telecommunication network provides count data. The estimation is assisted
#' with official data at a larger time scale (either from a population register or from a survey).
#'
#' The model rests on two working assumptions:
#' \itemize{
#'  \item Given that mobile phone data and official data operate at different time scales, we assume
#'  that there exists an initial time instant in which we can equate population figures from both
#'  sources.
#'
#'  \item The mobility patterns of individuals do not depend on the mobile network operator which
#'  they are subscribed to.
#'
#' }
#'
#' The model works in two stages. Firstly at the initial time instant, we use data from both sources
#' to make the inference for the actual population counts in each cell. Secondly, the time evolution
#' of these counts are produced using the transition matrices from cell to cell of individuals
#' provided by the mobile network operator.
#'
#' The essence of the model is to emulate the ecological sampling setting in which the number of
#' detected individuals in each cell follows a binomial distribution \eqn{Bin(N_{i}, p_{i})} whose
#' parameter \eqn{N_{i}} is the target of the model and is assigned a weakly informative prior and
#' the detection probability is also assigned a weakly informative prior based upon both data
#' sources.
#'
#' @section Computational paradigm:
#' Computations are conducted following the Bayesian paradigm. In this sense the generation of
#' simulated populations according to different probability distributions is at the core of the
#' package. In this sense the package contains basically three types of functions:
#'
#' \itemize{
#'
#'  \item Auxiliary functions, providing computation of mathematical functions such as the ratio of
#'  two beta functions, the confluent hypergeometric function, an optimization routine for a
#'  concrete probability distribution, etc. Examples of these functions are \code{\link{ratioBeta}},
#'   \code{\link{kummer}}, \code{\link{Phi}}, \code{\link{modeLambda}}.
#'
#'
#'  \item Distribution-relation functions, providing computation regarding the generation of random
#'  deviates according to different probability distributions comprising both priors, posteriors,
#'  and the generation of parameter specifications for these distributions. Examples of these
#'  functions are \code{\link{dtriang}}, \code{\link{rtriang}}, \code{\link{ptriang}},
#'  \code{\link{qtriang}}, \code{\link{dlambda}}, \code{\link{rlambda}}, \code{\link{rmatProb}},
#'  \code{\link{rN0}}, \code{\link{rNt}}, \code{\link{rNtcondN0}}, \code{\link{rg}},
#'  \code{\link{rp}}, \code{\link{alphaPrior}}, \code{\link{genAlpha}}, \code{\link{genUV}}.
#'
#'  \item Estimation-relation functions, providing computation of estimates based upon the
#'  populations generated with the preceding functions. Examples of these functions are
#'  \code{\link{postN0}}, \code{\link{postNt}}, \code{\link{postNtcondN0}}.
#'
#' }
#'
#' @docType package
#' @name pestim
NULL
