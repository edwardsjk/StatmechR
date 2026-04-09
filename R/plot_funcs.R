## PLOT FXNS ##

#' Plot full epidemic curves
#'
#' @description
#' `plot_true_curves()` plots the full epidemic curves from the data. Options to include or exclude a legend are available.
#'
#' @param dat A data frame providing the epidemic data over time for each location
#' @param X A character object providing the column name for the x-axis variable. This should be a variable that is in a Date format.
#' @param plot.group A character objecting providing the colunm name for the variable that the data will be grouped by
#' @param count A character object providing the column name for the count variable
#' @param legend A logical object indicating whether a legend key should be shown. The default is `TRUE`
#'
#' @return A ggplot of the epidemic curves over the x-axis variable and grouped by the plot_group variable
#'
#' @importFrom stringr str_to_title
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr ungroup
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
#'  plot_true_curves(dat = epidemic_data,
#'                   X = "epiweek_date",
#'                   plot_group = "district",
#'                   count = "cases", legend = TRUE)
#'
plot_true_curves <- function(dat,
                             X,
                             plot.group,
                             count,
                             legend = FALSE){

  rlang::check_required(dat)
  rlang::check_required(X)
  rlang::check_required(plot.group)
  rlang::check_required(count)

  check_empty(dat, "data.frame")
  check_class(dat, "data.frame")

  X_expr          <- rlang::enexpr(X)
  plot.group_expr <- rlang::enexpr(plot.group)
  count_expr      <- rlang::enexpr(count)

  check_single_col(X_expr,          arg = "X")
  check_single_col(plot.group_expr, arg = "plot.group")
  check_single_col(count_expr,      arg = "count")

  X_name <- rlang::as_name(rlang::ensym(X))
  plot.group_name <- rlang::as_name(rlang::ensym(plot.group))
  count_name <- rlang::as_name(rlang::ensym(count))

  check_contains_cols(X_name, dat, "X")
  check_contains_cols(plot.group_name, dat, "plot.group")
  check_contains_cols(count_name, dat, "count")

  check_NAs(dat[[X_name]], arg = "X", threshold = "all")
  check_NAs(dat[[count_name]], arg = "count", threshold = "all")

  check_class(dat[[count_name]], "numeric", "count")

  plotdat <- dat |>
    filter(!is.na({{count}})) |>
    group_by(pick(c({{plot.group}}, {{X}}))) |>
    summarize(cases = sum({{count}})) |>
    ungroup()

  if(legend == TRUE){

    tplot <- ggplot(plotdat, aes(x = {{X}}, y = cases,
                                 color = as.factor({{plot.group}}))) +
      geom_line() +
      labs(x = "Date", y = "New cases",
           color = rlang::as_name(rlang::ensym(plot.group))) +
      theme(axis.title.x = element_text(size = 20),
             axis.title.y = element_text(size = 20),
             axis.text.x = element_text(size = 18),
             axis.text.y = element_text(size = 18),
             legend.title = element_text(size = 20),
             legend.text = element_text(size = 18))

  }else{

    tplot <- ggplot(data = plotdat,
                    aes(x = {{X}},y = cases,
                        color = as.factor({{plot.group}}))) +
      geom_line() +
      labs(x = "Date", y = "New cases") +
      theme(axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.position = "none")

  }

  return(tplot)

}

#' Plot masked epidemic curves
#'
#' @description
#' `plot_mask_curves()` plots the masked curves of the epidemics, showing the curve up to the mask date in bold and the curve after the mask date in a faded/muted color. Options to include or exclude a legend are available.
#'
#' @param dat A data frame providing the full epidemic data over time for each location.
#' @param maskdat A data frame providing the masked epidemic data over time for each location.
#' @param X A character object providing the column name for the x axis variable. This variable should be in Date format
#' @param plot.group A character objecting providing the column name for the variable that the data will be grouped by.
#' @param count A character object providing the column name for the count variable
#' @param legend A logical object indicating whether a legend key should be shown. The default is `TRUE`.
#'
#' @return A ggplot of the masked epidemic curves over the x-axis variable and grouped by the plot_group variable.
#'
#' @importFrom stringr str_to_title
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
#' plot_mask_curves(dat = epidemic_data,
#'                  maskdat = masked_data,
#'                  X = "epiweek_date",
#'                  plot_group = "district",
#'                  count = "cases",
#'                  legend = TRUE)
#'
plot_mask_curves <- function(dat,
                             maskdat,
                             X,
                             plot.group,
                             count,
                             legend = FALSE){

  rlang::check_required(dat)
  rlang::check_required(maskdat)
  rlang::check_required(X)
  rlang::check_required(plot.group)
  rlang::check_required(count)

  check_empty(dat, "data.frame")
  check_empty(maskdat, "data.frame")
  check_class(dat, "data.frame")
  check_class(maskdat, "data.frame")

  X_expr          <- rlang::enexpr(X)
  plot.group_expr <- rlang::enexpr(plot.group)
  count_expr      <- rlang::enexpr(count)

  check_single_col(X_expr,          arg = "X")
  check_single_col(plot.group_expr, arg = "plot.group")
  check_single_col(count_expr,      arg = "count")

  X_name <- rlang::as_name(rlang::ensym(X))
  plot.group_name <- rlang::as_name(rlang::ensym(plot.group))
  count_name <- rlang::as_name(rlang::ensym(count))

  check_contains_cols(X_name, dat, "X")
  check_contains_cols(X_name, maskdat, "X")
  check_contains_cols(plot.group_name, dat, "plot.group")
  check_contains_cols(plot.group_name, maskdat, "plot.group")
  check_contains_cols(count_name, dat, "count")
  check_contains_cols(count_name, maskdat, "count")

  check_NAs(dat[[X_name]], arg = "X", threshold = "all")
  check_NAs(dat[[count_name]], arg = "count", threshold = "all")
  check_NAs(maskdat[[X_name]], arg = "X", threshold = "all")
  check_NAs(maskdat[[count_name]], arg = "count", threshold = "all")

  check_class(dat[[count_name]], "numeric", "count")
  check_class(maskdat[[count_name]], "numeric", "count")

  plotdat <- maskdat |>
    filter(!is.na({{count}})) |>
    group_by(pick({{plot.group}}, {{X}})) |>
    summarize(cases = sum({{count}})) |>
    ungroup() |>
    group_by(pick({{plot.group}})) |>
    filter(cases > 0) |>
    ungroup()

  plotdat_all <- dat |>
    filter(!is.na({{count}})) |>
    group_by(pick({{plot.group}}, {{X}})) |>
    summarize(cases = sum({{count}})) |>
    ungroup() |>
    group_by(pick({{plot.group}})) |>
    filter(cases > 0) |>
    ungroup()

  if(legend == TRUE){

    tplot <- ggplot() +
      geom_line(aes(x = {{X}},
                    y = cases, color = as.factor({{plot.group}})),
                data = plotdat) +
      geom_line(aes(x = {{X}}, y = cases,
                    color = as.factor({{plot.group}})),
                data = plotdat_all, alpha = 0.2) +
      labs(x = "Date", y = "New cases",
           lty = rlang::as_name(rlang::ensym(plot.group))) +
      theme(axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18))

  }else{

    tplot <- ggplot() +
      geom_line(aes(x = {{X}}, y = cases,
                    color = as.factor({{plot.group}})),
                data = plotdat) +
      geom_line(aes(x = {{X}}, y = cases,
                    color = as.factor({{plot.group}})),
                data = plotdat_all, alpha = 0.2) +
      labs(x = "Date", y = "New cases") +
      theme(axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.position = "none")


  }

  return(tplot)

}

#' Plot bias curves for the combined model
#'
#' @description
#' `plot_bias()` plots the error between the combined model's prediction and the true epidemic total at each iteration of the combined model. This allows the user to view changes in error and combined model performance over time.
#'
#' @param trueK A dataframe providing the true final epidemic sizes for each location
#' @param iter.results A dataframe providing the combined model final size estimate for each location at each iteration of the combined model.
#'
#' @return A ggplot object of the bias curves over the iterations of the combined model
#'
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
#' @importFrom dplyr arrange
#' @importFrom dplyr row_number
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
#' plot_bias(trueK = true_final,
#'           iter_results = Khist)
#'
plot_bias <- function(trueK,
                      iter.results){

  locs <- trueK |>
    ungroup() |>
    mutate(loc = row_number())

  bias_df <- data.frame(iter_results) |>
    mutate(iter = as.numeric(row.names(.))) |>
    pivot_longer(cols = -iter, names_to = "loc", values_to = "est",
                 names_prefix = "X") |>
    mutate(loc = as.numeric(gsub("[^[:digit:]]", "", loc))) |>
    left_join(locs) |>
    mutate(bias = est - K) |>
    arrange(loc, iter)

  bias_plot <- ggplot() +
    geom_line(aes(x = iter, y = bias, group = loc, color = factor(loc)),
              data = bias_df, alpha = 0.5) +
    labs(x = "Iteration", y = "Bias") +
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.position = "none")

  return(bias_plot)

}



