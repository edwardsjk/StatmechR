## UTILS FXNS ##

#' Extract the final epidemic size from the epidemic data
#'
#' @param dat A data frame providing the epidemic cases over time for each location
#' @param groups A character object providing the variable for the location level at which the final epidemic size should be calculated
#'
#' @return A data frame providing the location and the corresponding final epidemic size
#'
#' @importFrom dplyr ungroup
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#'
#' @export
#'
#' @examples
#'
#'  get_trueK(dat = epidemic_data, groups = "district")
#'
get_trueK <- function(dat, case_col, groups){

  grouping <- lapply(groups, function(x){

    which(names(dat) == x)

  })

  names(dat)[which(names(dat) == case_col)] <- "n_cases"

  K <- dat |>
    ungroup() |>
    group_by(do.call(pick, grouping)) |>
    summarise(K = sum(n_cases)) |>
    ungroup()

  return(K)

}

