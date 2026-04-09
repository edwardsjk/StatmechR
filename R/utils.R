## UTILS FXNS ##

#' Extract final epidemic size from the epidemic data
#'
#' @description
#' `get_trueK()` calculates the final epidemic size from epidemic data, given the column name containing the cases and the grouping variable.
#'
#'
#' @param dat A data frame providing the epidemic cases over time for each location
#' @param case.col A character object providing the variable for the cases
#' @param groups A character object providing the variable for the location level at which the final epidemic size should be calculated
#'
#' @return A data frame providing the location and the corresponding final epidemic size
#'
#' @importFrom dplyr ungroup
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom cli cli_abort
#'
#' @export
#'
#' @examples
#'
#'  get_trueK(dat = epidemic_data, case.col = n_cases,
#'            groups = district)
#'
get_trueK <- function(dat,
                      case.col,
                      groups){

  ## Checks ##

  rlang::check_required(dat)
  rlang::check_required(case.col)
  rlang::check_required(groups)

  check_empty(dat, "data.frame")
  check_class(dat, "data.frame")

  case_name  <- rlang::as_name(rlang::ensym(case.col))

  groups_expr <- rlang::enexpr(groups)
  group_names <- if (rlang::is_call(groups_expr)) {
    sapply(rlang::call_args(groups_expr), rlang::as_name)
  } else {
    rlang::as_name(groups_expr)
  }

  check_contains_cols(case_name,  dat, "case.col")
  check_contains_cols(group_names, dat, "groups")
  check_contains_cols(case_name, dat, "case.col")
  check_contains_cols(group_names, dat, "groups")
  check_class(dat[[case_name]], "numeric", "case.col")
  for (g in group_names) check_NAs(dat[[g]], arg = g,
                                   threshold = "any")

  ## Function ##

  K <- dat |>
    ungroup() |>
    dplyr::filter(!is.na({{case.col}})) |>
    group_by(pick({{ groups }})) |>
    summarise(K = sum({{case.col}})) |>
    ungroup()

  return(K)

}

