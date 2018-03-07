#' Dataset with simulated data for population counts.
#'
#' This dataset provides population counts moving from each pair of cells at succesive time instants
#'  for a simulated true population, a corresponding official population in a register and a
#'  population detected with a mobile telecommunication network.
#'
#' @format A \linkS4class{data.table} with 96768 rows and 6 variables:
#' \describe{
#'   \item{ID_CELL_INI}{identification code for each initial cell in the displacements}
#'   \item{ID_CELL_END}{identification code for each final cell in the displacements}
#'   \item{ID_T}{identification code of each time moment. It is very important to underline that the
#'   table collects always displacements between the initial time instant and the corresponding time
#'   instant specified by \code{ID_T}}
#'   \item{N_REG}{counts according to the population register. Note that these counts do not evolve
#'   in time}
#'   \item{N_0}{counts of the simulated true population}
#'   \item{N_MNO_1}{counts of individuals detected by the Mobile Network Operator}
#' }
"MobPop"
