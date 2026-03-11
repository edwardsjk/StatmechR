## PLOT FXNS ##

#' Plot the full epidemic curve
#'
#' @param dat A data frame providing the epidemic data over time for each location
#' @param X A character object providing the column name for the x-axis variable. This should be a variable that is in a Date format.
#' @param plot_group A character objecting providing the colunm name for the variable that the data will be grouped by
#' @param count A character object providing the column name for the count variable
#' @param legend A logical object indicating whether a legend key should be shown. The default is `TRUE`
#'
#' @return A ggplot of the epidemic curves over the x-axis variable and grouped by the plot_group variable
#'
#' @importFrom stringr str_to_title
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
#' @examples
#'
#'  plot_true_curves(dat = epidemic_data, X = "epiweek_date",
#'                   plot_group = "district", count = "cases", legend = TRUE)
#'
plot_true_curves <- function(dat, X, plot_group, count, legend = FALSE){

  names(dat)[which(names(dat) == plot_group)] <- "grouping"

  names(dat)[which(names(dat) == X)] <- "xaxis"

  names(dat)[which(names(dat) == count)] <- "case_count"

  plotdat <- dat |>
    group_by(grouping, xaxis) |>
    summarize(cases = sum(case_count)) |>
    ungroup()

  if(legend == TRUE){

    tplot <- ggplot(plotdat, aes(x = as.Date(x_axis), y = cases,
                                 color = as.factor(grouping))) +
      geom_line() +
      ylab("New cases") +
      xlab("Date") +
      labs(color = str_to_title(plot_group)) +
      gtheme(axis.title.x = element_text(size = 20),
             axis.title.y = element_text(size = 20),
             axis.text.x = element_text(size = 18),
             axis.text.y = element_text(size = 18),
             legend.title = element_text(size = 20),
             legend.text = element_text(size = 18))

  }else{

    tplot <- ggplot(data = plotdat,
                    aes(x = as.Date(xaxis),y = cases,
                        color = as.factor(grouping))) +
      geom_line() +
      ylab("New cases") +
      xlab("Date") +
      theme(axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.position = "none")

  }

  return(tplot)

}

#' Plot the masked epidemic curves
#'
#' @param dat A data frame providing the full epidemic data over time for each location.
#' @param maskdat A data frame providing the masked epidemic data over time for each location.
#' @param X A character object providing the column name for the xaxis variable. This variable should be in Date format
#' @param plot_group A character objecting providing the column name for the variable that the data will be grouped by.
#' @param count A character object providing the column name for the count variable
#' @param legend A logical object indicating whether a legend key should be shown. The default is `TRUE`.
#'
#' @return A ggplot of the masked epidemic curves over the x-axis variable and grouped by the plot_group variable.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stringr str_to_title
#'
#' @export
#'
#' @examples
#'
#' plot_true_curves(dat = epidemic_data, maskdat = masked_data,
#'                  X = "epiweek_date", plot_group = "district", count = "cases",
#'                  legend = TRUE)
#'
plot_mask_curves <- function(dat, maskdat, X, plot_group, count, legend = FALSE){

  names(dat)[which(names(dat) == plot_group)] <- "grouping"
  names(dat)[which(names(dat) == X)] <- "xaxis"
  names(dat)[which(names(dat) == count)] <- "case_count"

  names(maskdat)[which(names(maskdat) == plot_group)] <- "grouping"
  names(maskdat)[which(names(maskdat) == X)] <- "xaxis"
  names(maskdat)[which(names(maskdat) == count)] <- "case_count"

  plotdat <- mask_dat |>
    group_by(grouping, xaxis) |>
    summarize(cases = sum(case_count)) |>
    ungroup() |>
    group_by(grouping) |>
    filter(cases > 0) |>
    ungroup()

  plotdat_all <- dat |>
    group_by(grouping, xaxis) |>
    summarize(cases = sum(case_count)) |>
    ungroup() |>
    group_by(grouping) |>
    filter(cases > 0) |>
    ungroup()

  if(legend == TRUE){

    tplot <- ggplot() +
      geom_line(aes(x = as.Date(xaxis),
                    y = cases, color = as.factor(grouping)),
                data = plotdat) +
      geom_line(aes(x = as.Date(xaxis), y = cases,
                    color = as.factor(grouping2)),
                data = plotdat_all, alpha = 0.2) +
      ylab("New cases") +
      xlab("Date") +
      labs(lty = str_to_title(plot_group)) +
      theme(axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18))

  }else{

    tplot <- ggplot() +
      geom_line(aes(x = as.Date(xaxis), y = cases,
                    color = as.factor(grouping1)),
                data = plotdat) +
      geom_line(aes(x = as.Date(xaxis), y = cases,
                    color = as.factor(grouping)),
                data = plotdat_all, alpha = 0.2) +
      ylab("New cases") +
      xlab("Date") +
      theme(axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.position = "none")


  }

  return(tplot)

}

#' Plot the bias curves for the combined model
#'
#' @param trueK A dataframe providing the true final epidemic sizes for each location
#' @param iter_results A dataframe providing the combined model final size estimate for each location at each iteration of the combined model.
#'
#' @return A ggplot object of the bias curves over the iterations of the combined model
#'
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
#' plot_bias(trueK = true_final, iter_results = Khist)
#'
plot_bias <- function(trueK, iter_results){

  locs <- trueK %>%
    ungroup() %>%
    mutate(loc = row_number())

  bias_df <- data.frame(iter_results) %>%
    mutate(iter = as.numeric(row.names(.))) %>%
    pivot_longer(cols = -iter, names_to = "loc", values_to = "est",
                 names_prefix = "X") %>%
    mutate(loc = as.numeric(gsub("[^[:digit:]]", "", loc))) %>%
    left_join(locs) %>%
    mutate(bias = est - K)%>%
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



