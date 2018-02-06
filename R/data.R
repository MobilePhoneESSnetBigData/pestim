#' True population in each cell at each time step.
#'
#' A dataset containing the true population for a set of cells for succesive moments of time.
#'
#' @format A data table with 60 rows and 3 variables:
#' \describe{
#'   \item{Time}{\code{int}, time steps for which we want to estimate the population in each cell}
#'   \item{Cell}{\code{factor}, the cell number for which we estimate the population}
#'   \item{N0True}{\code{numeric}, the true population in each cell, at each time moment}
#' }
"TruePop"

#' Simulated population mobility.
#'
#' A dataset containing simulated population mobility from one cell to another cell.
#'
#' @format A data table with 760 rows and 6 variables:
#' \describe{
#'   \item{ID_CELL_INI}{\code{int}, time steps for which we want to estimate the population in each cell}
#'   \item{ID_CELL_FIN}{\code{int}, time steps for which we want to estimate the population in each cell}
#'   \item{ID_T}{\code{int}, time steps for which we want to estimate the population in each cell}
#'   \item{N_REG}{\code{numeric}, the number of population in each cell from the population administrative register}
#'   \item{N0}{\code{numeric}, the initial number of population in each cell}
#'   \item{N_MNO_1}{\code{numeric}, the number of population as given by the mobile phone operator}
#' }
"Sim.PopMob"
