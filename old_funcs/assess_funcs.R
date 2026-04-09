## ASSESSMENT FXNS ##

#' Calculate cases averted by allocation of a vaccine
#'
#' @description
#' `get_cases_prevented()` calculates the cases that could be averted if allocation of a vaccine was prioritized using the total cases expected. Specifically it uses the total number of cases expected and the current number of cases to determine which location have the most cases remaining.
#'
#'
#' @param method.pred A mathematical expression used to estimate the number of remaining cases
#' @param method.name A character object providing the name of the method
#' @param dat A dataframe providing the location name (`loc`), the true epidemic size (`K`), the observed cases (`observed`), the predicted sizes from the models, and the population size for the location (`N`).
#' @param maxvax A numeric object providing the max number of vaccine doses that can be allocated. The default is `3e6`.
#' @param coverage A numeric object providing the coverage of the vaccine as a proportion (i.e. between 0 and 1).
#' @param vax.eff A numeric object providing the vaccine efficacy as a proportion (between 0 and 1).
#'
#' @returns A dataframe containing the following columns, that shows where
#'          doses were allocated and the number of infections prevented:
#'            - `loc` which provides the location name.
#'            - `method` which provides the name of the method used to
#'              calculated the estimated remaining cases.
#'            - `cumv2` which provides the cumulative number of vaccine doses
#'              allocated after each location.
#'            - `inf_prevented` which provides the cumulative number of
#'              infections/cases averted
#'
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr select
#'
#' @export
#'
#' @examples
#'
#' get_cases_prevented(method.pred = stat - observed, method.name = "stat",
#'                     dat = my_dataset, maxvax = 3e6, coverage = 0.75,
#'                     vax_eff = 0.75)
#'
get_cases_prevented <- function(method.pred, method.name,
                                dat = flat_res, maxvax = 3e6,
                                coverage, vax.eff) {

  res <- dat %>%
    ungroup() %>%
    mutate(method = method.name,
           u = runif(n = nrow(dat)),
           est_remaining = {{method.pred}},
           true_remaining = K - observed,
           sortval = ifelse(method == "random", u, (est_remaining)/N)) |>
    arrange(-sortval) |>
    mutate(vaxxed = (N) * coverage,
           cumv = cumsum(vaxxed)) |>
    filter(lag(cumv, n = 1, default = 0) <= maxvax) |>
    mutate(cumv2 = pmin(cumv, maxvax),
           vaxhere = ifelse(cumv2 == cumv, vaxxed, cumv - cumv2),
           inf_prevented = cumsum((true_remaining * vaxhere / vaxxed) * coverage * vax.eff)) |>
    select(loc, method, cumv2, inf_prevented)

  return(res)

}
