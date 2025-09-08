#### NECESSARY FUNCTIONS ####

## Required packages: tibble; tidyverse

###################################
#### Data processing functions ####
###################################

#' Simulate artificial covariates for locatons
#'
#' @param N A numeric object providing the number of locations that should be simulated.
#' @param poverty A character object providing the level of poverty in the area. Input options are `"high"`, `"average"` or `"low"`. The default is `"average"`.
#' @param climate A numeric vector between 0 and 1 providing the proportion of locations that should be assigned to the following three climate types: "wet", "moderate" and "arid". The default is 0.33 for each climate type
#' @param temperature A numeric vector providing temperature (in celcius) associated with the three different climate types. The default is `c(30, 25, 25)`.
#' @param water_san A character object providing the level of water and sanitation in the area. Input options are `"high"`, `"average"` or `"low"`.  The default is `"average"`.
#' @param stand_water A character object providing the level of standing water in the area.
#' Input options are `"high"`, `"average"` or `"low"`.  The default is `"average"`.
#' @param education A character object providing the level of education in the area. Input options are `"high"`, `"average"` or `"low"`.  The default is `"average"`.
#'
#' @return A dataframe with a row for each location and a column for each covariate.
#'
#' @importFrom dplyr tibble
#' @importFrom dplyr mutate
#'
#' @export
#'
#' @examples
#' sim_covariates(N = 300,
#'                poverty = "high",
#'                climate = c(0.3, 0.4, 0.4),
#'                temperature = c(40, 35, 50),
#'                water_san = "low",
#'                stand_water = "average",
#'                education = "low")
#'
sim_covariates <- function(N,
                            poverty = "average",
                            climate = c(0.33, 0.33, 0.33),
                            temperature = c(30, 20, 35),
                            water_san = "average",
                            stand_water = "average",
                            education = "average") {

  watsan_loc <- ifelse(water_san == "average", 0, ifelse(water_san == "low", 1.5, -1.5))

  stand_loc <- ifelse(stand_water == "average", 0, ifelse(stand_water == "low", 1.5, -1.5))

  edu_loc <- ifelse(education == "average", 0, ifelse(education == "low", 1.5, -1.5))

  if(poverty == "average"){

    rc <- tibble(pov = rbeta(N, 3, 3, 0.5))

  }

  if(poverty == "low"){

    rc <- tibble(pov = rbeta(N, 1, 3, 0))

  }

  if(poverty == "high"){

    rc <- tibble(pov = rbeta(N, 3, 1, 0))

  }

  rc <- rc |> mutate(
    clim = sample(c("wet","mod","arid"), N, replace = T, prob = climate), #climate
    temp = rnorm(N, (clim == "wet") * temperature[1] + (clim == "mod") * temperature[2] + (clim == "arid") * temperature[3], 5), #Temprature
    watsan = plogis(rnorm(N, 2 + -2 * pov), location = watsan_loc, scale = 0.5), #wat/san access
    stand = plogis(rnorm(N, -3 + 2 * (clim == "wet")+ -2 * (clim == "dry") + 2 * (clim == "wet" | clim == "mod")), location = stand_loc, scale = 0.5),#standing water
    educ = plogis(rnorm(N, 1 + -2 * pov), location = edu_loc, scale = 0.5) #Education
  )

  return(rc)

}

##

#' Compute R0 values from simulated covariate data
#'
#' @param covs Dataframe providing the covariate values for each location
#' @param b0 Numeric object providing the R0 intercept on the log scale. The default is `log(2)`
#' @param sd Numeric object providing the standard deviation in R0 on the log scale. The default is `0.1`.
#' @param b1 Numeric object providing the coefficient for the wat/san access covariate. The default is `log(0.75)`.
#' @param b2 Numeric object providing the coefficient for the standing water covariate. The default is `log(1.2)`
#' @param b3 Numeric object providing the coefficient for the interaction of wat/san access and standing water covariates. The default is `log(1.2)`
#' @param temp_scale Numeric object providing a scale parameter for a logistic distribution of the temperature covariate. The default is `4`
#' @param temp_loc Numeric object providing a location parameter for a logistic distribution of the temperature covariate. The default is `20`
#' @param b4 Numeric object providing the coefficient for the temperature covariate. The default is `log(1.2)`
#'
#' @return A numeric vector providing the calculated R0 values for each location
#'
#' @importFrom dplyr mutate
#'
#' @export
#'
#' @examples
#'
#' compute_Rs(covariate_data, b0 = log(2), sd = .1,
#'            b1 = log(.75), b2 = log(1.2), b3 = log(1.2),
#'            temp_scale = 4, temp_loc = 20, b4 = log(1.2))
#'
compute_Rs <- function(covs, b0 = log(2), sd = .1, b1 = log(.75),
                       b2 = log(1.2), b3=log(1.2),
                       temp_scale = 4, temp_loc=20,
                       b4 = log(1.2) ) {

  ##first normalize the covariates
  covs <- covs |> mutate(watsan_stand = ((1-watsan)*stand-mean((1-watsan)*stand))/sd((1-watsan)*stand), #normalise interactoin term
                         watsan = (watsan-mean(watsan))/sd(watsan),
                         stand = (stand-mean(stand)/sd(stand)),
                         R0 = exp(rnorm(n(),b0+b1*watsan+b2*stand+b3*watsan_stand +
                                          b4*plogis(temp, temp_loc, temp_scale), sd)))

  return(covs$R0)

}

##

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
get_trueK <- function(dat, groups){

  grouping <- lapply(groups, function(x){

    which(names(dat) == x)

  })

  K <- dat |>
    ungroup() |>
    group_by(do.call(pick, grouping)) |>
    summarise(K = sum(n_cases)) |>
    ungroup()

  return(K)

}

##

#' Plot the full epidemic curve
#'
#' @param dat A data frame providing the epidemic data over time for each location
#' @param X A character object providing the variable that should be on the x-axis
#' @param plot_group A character objecting providing the grouping variable for the data
#' @param legend A logical object indicating whether a legend key should be shown. The default is `TRUE`
#'
#' @return A ggplot of the epidemic curves over the x-axis variable and grouped by the plot_group variable
#'
#' @importFrom stringr str_to_title
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr ungroup
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
#'  plot_true_curves(dat = epidemic_data, X = "epiweek_date",
#'                   plot_group = "district", legend = TRUE)
#'
plot_true_curves <- function(dat, X, plot_group, legend = FALSE){

  grouping <- unlist(lapply(plot_group, function(b){

    which(names(dat) == b)

  }))

  xaxis <- unlist(lapply(X, function(c){

    which(names(dat) == c)

  }))

  plotdat <- dat |>
    group_by(.[[grouping]], .[[xaxis]]) |>
    summarize(cases = sum(n_cases)) |>
    ungroup()

  if(legend == TRUE){

    tplot <- ggplot(plotdat, aes(x = as.Date(`.[[xaxis]]`), y = cases,
                                          color = as.factor(`.[[grouping]]`))) +
      geom_line() +
      ylab("New cases") +
      xlab("Week") +
      labs(color = str_to_title(plot_group)) +
      gtheme(axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18))

  }else{

    tplot <- ggplot(data = plotdat,
                    aes(x = as.Date(`.[[xaxis]]`),y = cases,
                        color = as.factor(`.[[grouping]]`))) +
      geom_line() +
      ylab("New cases") +
      xlab("Week") +
      theme(axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.position = "none")

  }

  return(tplot)

}

##

#' Plot the masked epidemic curves
#'
#' @param dat A data frame providing the full epidemic data over time for each location.
#' @param maskdat A data frame providing the masked epidemic data over time for each location.
#' @param X A character object providing the variable that should be on the x-axis.
#' @param plot_group A character objecting providing the grouping variable for the data.
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
#'                  X = "epiweek_date", plot_group = "district", legend = TRUE)
#'
plot_mask_curves <- function(dat, maskdat, X, plot_group, legend = FALSE){

  grouping1 <- unlist(lapply(plot_group, function(x){

    which(names(maskdat) == x)

  }))

  grouping2 <- unlist(lapply(plot_group, function(x){

    which(names(dat) == x)

  }))

  xaxis <- unlist(lapply(X, function(c){

    which(names(dat) == c)

  }))

  plotdat <- mask_dat |>
    group_by(.[[grouping1]], .[[xaxis]]) |>
    summarize(cases = sum(n_cases)) |>
    ungroup() |>
    group_by(.[[grouping1]]) |>
    filter(cases > 0) |>
    ungroup()

  plotdat_all <- dat |>
    group_by(.[[grouping2]], .[[xaxis]]) |>
    summarize(cases = sum(n_cases)) |>
    ungroup() |>
    group_by(.[[grouping2]]) |>
    filter(cases > 0) |>
    ungroup()

  if(legend == TRUE){

    tplot <- ggplot() +
      geom_line(aes(x = as.Date(`.[[xaxis]]`),
                    y = cases, color = as.factor(`.[[grouping1]]`)),
                data = plotdat) +
      geom_line(aes(x = as.Date(`.[[xaxis]]`), y = cases,
                                      color = as.factor(`.[[grouping2]]`)),
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
      geom_line(aes(x = as.Date(`.[[xaxis]]`), y = cases,
                    color = as.factor(`.[[grouping1]]`)),
                data = plotdat) +
      geom_line(aes(x = as.Date(`.[[xaxis]]`), y = cases,
                    color = as.factor(`.[[grouping2]]`)),
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

##

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

############################
#### Modeling functions ####
############################

#' Simulate a Gaussian epidemic curve
#'
#' @param pars A vector providing, in order:
#' - the estimated epidemic size
#' - the estimated peak time of the epidemic
#' - the estimated spread of the epidemic
#' @param times A numeric object providing the number of timesteps that the epidemic should be simulated for
#'
#' @return A vector of incidence values for each timestep
#' @export
#'
#' @examples
#'
#' normmdl(pars = c(250, 15, 5), time = 30)
#'
normmdl <- function(pars, times) {

  fs_i <- unname(unlist(pars[1]))

  peaktime_i <- unname(unlist(pars[2]))

  spread_i <- exp(unname(unlist(pars[3])))

  normpred <- fs_i * (pnorm(1:times, peaktime_i, spread_i) -
                        pnorm((1:times) - 1, peaktime_i, spread_i))

  return(normpred)
}

##

#' Train a SuperLearner ensemble in parallel
#'
#' @param x A dataframe providing the variables for each location.
#' @param y A numeric vector providing the outcome for each location.
#' @param cores A numeric object providing the number of cores to be used in the parallel process. The defaultt value is `1`.
#' @param family A character object providing the error distribution. Currenly only `gaussian` (the default) and `poisson` are supported.
#'
#' @return A `SuperLearner` object containing the trained model
#'
#' @importFrom SuperLearner snowSuperLearner
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterSetRNGStream
#' @importFrom parallel stopCluster
#'
#' @export
#'
#' @examples
#'
#' stat.mdl.sl.fit.para(x = my_vars, y = epi_size,
#'                      cores = 4, family = "gaussian)
#'
stat.mdl.sl.fit.para <- function(x, y, cores = 1, family = "gaussian") {

  cluster = makeCluster(cores)
  clusterEvalQ(cluster,
               {library(SuperLearner)
                 library(gam)
                 library(rpart)
                 library(randomForest)
                 library(e1071)          # For SL.svm
                 library(bartMachine)    # For SL.bartMachine
                 options(mc.cores = 1)})
  clusterSetRNGStream(cluster, 1)

  if(family == "gaussian"){

    rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "gaussian", cluster = cluster,
                           SL.library = list(c("SL.rpart", "screen.glmnet"),
                                             c("SL.randomForest", "screen.glmnet"),
                                             c("SL.glm", "screen.glmnet"),
                                             c('SL.gam', "screen.glmnet"),
                                             c("SL.glmnet", "screen.glmnet"),
                                             c("SL.xgboost", "screen.glmnet"),
                                             c("SL.svm", "screen.glmnet"),
                                             "SL.glmnet"))

    stopCluster(cluster)
    gc()

  }

  if(family == "poisson"){

    rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "poisson",
                           cluster = cluster,
                           SL.library = list(c("SL.glm", "screen.glmnet"),
                                             c('SL.gam', "screen.glmnet"),
                                             c("SL.glmnet", "screen.glmnet"),
                                             "SL.glmnet",
                                             c("SL.bartMachine", "screen.glmnet")))

    stopCluster(cluster)
    gc()
  }

  return(rc)

}

#' Train a SuperLearner ensemble in parallel
#'
#' @param x A dataframe providing the variables for each location.
#' @param y A numeric vector providing the outcome for each location.
#' @param cores A numeric object providing the number of cores to be used in the parallel process. The defaultt value is `1`.
#' @param family A character object providing the error distribution. Currenly only `gaussian` (the default) and `poisson` are supported.
#'
#' @return A `SuperLearner` object containing the trained model
#'
#' @importFrom SuperLearner snowSuperLearner
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterSetRNGStream
#' @importFrom parallel stopCluster
#'
#' @export
#'
#' @examples
#'
#' stat.mdl.sl.fit.para(x = my_vars, y = epi_size,
#'                      cores = 4, family = "gaussian)
#'
stat.mdl.sl.fit.para2 <- function(x, y, cores = 1, family = "gaussian") {

  cluster = makeCluster(cores)
  clusterEvalQ(cluster,
               {library(SuperLearner)
                library(gam)
                library(rpart)
                library(randomForest)
                library(e1071) # For SL.svm
                 options(mc.cores = 1)})
  clusterSetRNGStream(cluster, 1)

  if(family == "gaussian"){

    rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "gaussian", cluster = cluster,
                           SL.library = list(c("SL.rpart", "screen.glmnet"),
                                             c("SL.randomForest", "screen.glmnet"),
                                             c("SL.glm", "screen.glmnet"),
                                             c('SL.gam', "screen.glmnet"),
                                             c("SL.glmnet", "screen.glmnet"),
                                             c("SL.xgboost", "screen.glmnet"),
                                             c("SL.svm", "screen.glmnet"),
                                             "SL.glmnet"))

    stopCluster(cluster)
    gc()

  }

  if(family == "poisson"){

    rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "poisson",
                           cluster = cluster,
                           SL.library = list(c("SL.glm", "screen.glmnet"),
                                             c('SL.gam', "screen.glmnet"),
                                             c("SL.glmnet", "screen.glmnet"),
                                             "SL.glmnet"))

    stopCluster(cluster)
    gc()
  }

  return(rc)

}

##

#' Title
#'
#' @param x A dataframe providing the variables for each location.
#' @param y A numeric vector providing the outcome for each location.
#' @param family A character object providing the error distribution. Currenly only `gaussian` (the default) and `poisson` are supported.
#'
#' @return A `SuperLearner` object containing the trained model
#'
#' @importFrom SuperLearner SuperLearner
#'
#' @export
#'
#' @examples
#'
#' stat.mdl.sl.fit(x = my_vars, y = epi_size, family = "gaussian)
#'
stat.mdl.sl.fit <- function(x, y, family = "gaussian") {

  if(family == "gaussian"){

    rc <- SuperLearner(X = x, Y = y, family = "gaussian",
                       SL.library = list(c("SL.rpart", "screen.glmnet"),
                                         c("SL.randomForest", "screen.glmnet"),
                                         c("SL.glm", "screen.glmnet"),
                                         c('SL.gam', "screen.glmnet"),
                                         c("SL.glmnet", "screen.glmnet"),
                                         c("SL.xgboost", "screen.glmnet"),
                                         c("SL.svm", "screen.glmnet"),
                                         "SL.glmnet"))

  }

  if(family == "poisson"){

    rc <- SuperLearner(X = x, Y = y, family = "poisson",
                       SL.library = list(c("SL.glm", "screen.glmnet"),
                                         c('SL.gam', "screen.glmnet"),
                                         c("SL.glmnet", "screen.glmnet"),
                                         "SL.glmnet",
                                         c("SL.bartMachine", "screen.glmnet")))

  }

  return(rc)

}

##

#' Predict final epidemic sizes using a trained SuperLearner
#'
#' @param mdl A `SuperLearner` object providing the trained ensemble model.
#' @param x A dataframe providing the variables/features for each location.
#'
#' @return A numeric vector of predicted epidemic final sizes for each location.
#'
#' @export
#'
#' @examples
#'
#' stat.mdl.sl.pred(mdl = SL_model, x = my_vars)
#'
stat.mdl.sl.pred <- function(mdl, x) {

  pred <- pmax(0, predict(mdl, onlySL = T, newdata = x)$pred)

  return(pred)

}

##

#' Gaussian model wrapper function
#'
#' @param ecs A list object providing the observed or masked epidemic curves. Each element of the list object is a numeric vector.
#' @param pop_N A numeric vector providing the population size for each location
#' @param strt_vals A data frame providing the starting parameter values for the Gaussian model. These are the estimated final size, peak time, and spread. Parameters are represented by columns and the locations are represented by rows.
#' @param errorfxn A function object providing the error calculation for the optimization.
#' @param penaltyfunc A function object providing the penalty calculation used to penalize the optimization based `priorval`. The default is `NULL`.
#' @param priorval A numeric vector providing the prior values produced by the statistical model, which the optimization will be penalized towards. The default is `NULL`.
#' @param cores A numeric object providing the number of cores to be assigned to the model process when run in parallel.
#' @param tau A numeric object providing the number of time steps that the mechanistic model should be simulating.
#' @param timestep A numeric object providing the timestep frequency, in number of days, i.e. one week would be `7` while one day would be `1`.
#'
#' @return A list object containing the new parameter values produced by the optimization (for each location), and the convergence results of the optimization.
#'
#' @export
#'
#' @examples
#'
#' norm_em(ecs = my_curves, pop_N = population, errorfxn = norm_error,
#' penaltyfunc = NULL, priorval = NULL, cores = 1, tau = 52, timestep = 7)
#'
norm_em <- function(ecs, pop_N, strt_vals, errorfxn, penaltyfunc, priorval,
                    cores, tau, timestep) {

  model_fit <- fit_norm_model(ecs, normmdl, strt_vals, errorfxn,
                              priorfunc = penaltyfunc, prior = priorval,
                              cores, tau)

  return(model_fit)
}

##

#' Combined mechanistic and statistic iterative model
#'
#' @param epi_curves A list object providing the observed or masked epidemic curves. Each element of the list object is a numeric vector.
#' @param covdat A data frame object providing the covariate values for each location.
#' @param pop_N A numeric vector providing the population size for each location.
#' @param initK ## removing argument ##
#' @param epimdlfit A function object providing the wrapper used to run the chosen mechanistic model.
#' @param starting_vals A data frame providing the starting parameter values for the mechanistic model.
#' @param error_func A function object providing the error calculation for the optimization of the mechanistic model.
#' @param penalty_func A function object providing the penalty calculation used to penalize the mechanistic model optimization towards the results from the statistical model.
#' @param statmdlfit A function object providing the function used to fit the statistical model component.
#' @param statmdlpred A function object providing the function used to predict the outcome using the trained statistical model.
#' @param threshold A numeric object providing a threshold for the iteration difference. The iteration difference is difference between the current iteration prediction and previous iteration prediction. Once this drops below the threshold among all locations, the iterative process will stop. The default is `20`
#' @param max.iter A numeric object providing the maximum number of iterations that should be completed. If the threshold has not been met by this number of iterations, the iterative process will stop. The default is `100`
#' @param cores A numeric object providing the number of cores that should be assigned to run the function in parallel.
#' @param stat.family A character object providing the error distribution family for the statistical model. The options are `"gaussian"` and `"poisson"`.
#' @param tau A character object providing the number of timesteps that should be simulated for the mechanistic model prediction.
#' @param timestep A numeric object providng the number of days in each time step. For example, a weekly time step would be `timwestep = 7`
#'
#' @return A list object containing the following objects:
#' - `K` Most recent prediction of epidemic size from combined model
#' - `Kmech` Most recent prediction of epidemic size from the mechanistic component
#' - `Khist` Epidemic size predictions from the combined model at each iteration
#' - `Khist.mech` Epidemic size predictions from mechanistic component at each iteration
#' - `epi.params` Most recently optimized model parameters
#' - `params` Optimized model parameters for each iteration
#' - `converged` Convergence results from the most recent iteration
#' - `diff` Iteration difference between the last iteration and the one previous
#'
#'
#'
#' @export
#'
#' @examples
#'
#' em_func_model(epi_curves = my_curves, covdat = env_data, init_K = kmech,
#'               starting_vals = values, error_func = norm_error,
#'               statmdlfit = my_stat_fit, statmdlpred = my_stat_pred,
#'               threshold = 5, max.iter = 1000, cores = 20,
#'               stat.family = "gaussian", tau = 52, timestep = 7)
#'
em_func_model <- function(epi_curves, covdat, pop_N, initK, epimdlfit,
                          starting_vals, error_func, penalty_func, statmdlfit,
                          statmdlpred, threshold = 20, max.iter = 100, cores,
                          stat.family = "gaussian", tau = 52, timestep = NULL) {

  iter <- 0                   # initialize iter
  iter_diff <- 2 * threshold  # set iter_diff > threshold for first iter

  obs <- unname(unlist(lapply(epi_curves, function(x){sum(x)})))

  ##matrices for holding results.
  K <- matrix(nrow = max.iter, ncol = length(epi_curves))
  Kmech <- matrix(nrow = max.iter, ncol = length(epi_curves))

  prev_epi_mdl <- starting_vals

  ## Loop until ending criteria is met
  while (any(iter_diff >= threshold) & (iter < max.iter)) {
    print(iter)
    iter <- iter + 1                # update iter
    if(iter == 1){

      lastK <- initK

    }else{  # set lastK to starting values for first iter

      lastK <- K[iter - 1, ]

      } # otherwise, lastK is K from last iter

    ## Force lastK and first parameter to be >= obs

    lastK[which(lastK < obs)] <- obs[which(lastK < obs)]

    prev_epi_mdl[which(prev_epi_mdl[,1] < obs), 1] <- obs[which(prev_epi_mdl[,1] < obs)]

    if(iter == 1){

      fitepimdl <- epimdlfit(epi_curves, pop_N,
                             strt_vals = prev_epi_mdl,
                             error_func, penaltyfunc = NULL,
                             priorval = NULL, cores,
                             tau, timestep)
    }else{

      fitepimdl <- epimdlfit(epi_curves, pop_N,
                             strt_vals = prev_epi_mdl,
                             error_func, penalty_func,
                             priorval = lastK, cores,
                             tau, timestep)

    }

    params <- prev_epi_mdl

    tmp <- lapply(fitepimdl, function(x){x[[1]]})

    prev_epi_mdl <- as.data.frame(do.call(rbind, tmp))

    Kmech[iter, ] <- prev_epi_mdl[, 1]

    ## Fit the statistical model on this iteration

    Kmech2 <- Kmech[iter,]/(pop_N/1000)

    fitstatmdl <- statmdlfit(covdat, Kmech2, cores, stat.family)
    K[iter, ] <- statmdlpred(fitstatmdl, covdat) * (pop_N/1000)

    ## Check iter_diff
    if (iter > 1) {
      iter_diff <- abs(K[iter - 1, ] - K[iter, ])
    }

    print(rbind(lastK[1:5], Kmech[iter, 1:5], K[iter, 1:5]))

    gc()

  }

  converged <- lapply(fitepimdl, function(x){x[[2]]})

  return(list(K = K[iter, ],
              Kmech = Kmech[iter, ],
              Khist = K[1:iter, ],
              Khist.mech = Kmech[1:iter, ],
              epi.parms = fitepimdl,
              param.all = params,
              converged = converged,
              diff = iter_diff))

}

##

#' Optimization of a gaussian epidemic model
#'
#' @param ecs  A list object providing the observed or masked epidemic curves. Each element of the list object is a numeric vector.
#' @param epi_mdl_func A function object providing the function used to simulate the mechanistic model.
#' @param epi_mdl_pars A dataframe providing the parameter values for each parameter (columns) and location (rows).
#' @param error_func A function object providing the error calculation for the optimization of the mechanistic model.
#' @param priorfunc A function object providing the penalty calculation used to penalize the mechanistic model optimization towards the results from the statistical model.
#' @param prior A numeric vector providing the results from the statistical model that the `priorfunc` will penalize towards.
#' @param cores A numeric object providing the number of cores that should be assigned to run the function in parallel.
#' @param tau A character object providing the number of timesteps that should be simulated for the mechanistic model prediction.
#' @param timestep A numeric object providng the number of days in each time step. For example, a weekly time step would be `timestep = 7`
#'
#' @returns A list object consisting of, for each location:
#'            - a dataframe with the updated parameter values
#'            - a numeric object indicating whether the optimization
#'              converged. A `0` indicates convergence.
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#'
#' @export
#'
#' @examples
#'
#' fit_norm_model(ecs = my_curves, epi_mdl_func = my_epi_model,
#'                epi_mdl_pars = starting_values, error_func = norm_error,
#'                prior_func = NULL, prior = NULL, cores = 1, tau = 70,
#'                timestep = NULL)
#'
fit_norm_model <- function(ecs, epi_mdl_func, epi_mdl_pars,
                           error_func, priorfunc = NULL, prior = NULL,
                           cores, tau, timestep = NULL) {

  obj_fxn <- function(pars, ec, priorfunc = NULL, prior = NULL) {

    pred_curve <- epi_mdl_func(pars, length(ec))

    if(pars[1] == 0){

      pars[1] <- 1

    }

    if(!is.null(priorfunc)) {

      err <- error_func(ec, pred_curve, (pars[1])) +
        priorfunc((pars[1]), prior)

    }

    if(is.null(priorfunc)) {

      err <- error_func(ec, pred_curve, (pars[1]))

    }


    return(-err)

  }

  cluster <- makeCluster(cores)
  registerDoParallel(cluster)

  mod_res <- foreach(i = seq_along(ecs), .packages = c("dplyr", "StatmechR"))%dopar%{

    ec <- ecs[[i]]
    mdl_pars <- unname(unlist(epi_mdl_pars[i,]))

    if (!is.null(priorfunc)) {

      tmp <- optim(mdl_pars, obj_fxn, ec = ec, priorfunc = priorfunc,
                   prior = prior[i], method = "SANN", control = list(maxit = 1000))

    }

    if (is.null(priorfunc)) {

      tmp <- optim(mdl_pars, obj_fxn, ec = ec, method = "SANN", control = list(maxit = 1000))

    }

    converged <- tmp$convergence

    new_pars <- tmp$par

    return(list(new_pars, converged))

    gc()

  }

  stopCluster(cl = cluster)

  return(mod_res)
}

##

#' Optimization of a compartmental epidemic model
#'
#'
#' @param ecs  A list object providing the observed or masked epidemic curves. Each element of the list object is a numeric vector.
#' @param N A numeric vector providing the population size in each location.
#' @param epi_mdl_func A function object providing the function used to simulate the mechanistic model.
#' @param epi_mdl_pars A dataframe providing the parameter values for each parameter (columns) and location (rows).
#' @param error_func A function object providing the error calculation for the optimization of the mechanistic model.
#' @param priorfunc A function object providing the penalty calculation used to penalize the mechanistic model optimization towards the results from the statistical model.
#' @param prior A numeric vector providing the results from the statistical model that the `priorfunc` will penalize towards.
#' @param cores A numeric object providing the number of cores that should be assigned to run the function in parallel.
#' @param tau A character object providing the number of timesteps that should be simulated for the mechanistic model prediction.
#' @param timestep A numeric object providng the number of days in each time step. For example, a weekly time step would be `timestep = 7`
#'
#' @returns A list object consisting of, for each location:
#'            - a dataframe with the updated parameter values
#'            - a numeric object indicating whether the optimization
#'              converged. A `0` indicates convergence.
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#'
#' @export
#'
#' @examples
#'
#'fit_epi_model(ecs = my_curves, epi_mdl_func = my_epi_model,
#'              epi_mdl_pars = starting_values, error_func = norm_error,
#'              prior_func = NULL, prior = NULL, cores = 1, tau = 70,
#'              timestep = 7)
#'
fit_epi_model <- function(ecs, N, epi_mdl_func, epi_mdl_pars,
                          error_func, priorfunc = NULL, prior = NULL, cores,
                          tau, timestep) {

  obj_fxn <- function(pars, ec, N, priorfunc = NULL, prior = NULL, timestep) {

    pred_curve <- epi_mdl_func(N, pars[2:length(pars)], timestep, (length(ec) + 1))

    if(pars[1] == 0){

      pars[1] <- 1

    }

    if(!is.null(priorfunc)) {

      err <- error_func(ec, pred_curve$incident, pars[1]) +
        priorfunc((pars[1]), prior)
    }

    if(is.null(priorfunc)) {

      err <- error_func(ec, pred_curve$incident, pars[1])

    }

    return(-err)

  }

  cluster <- makeCluster(cores)
  registerDoParallel(cluster)

  mod_res <- foreach(i = seq_along(ecs), .packages = c("dplyr", "StatmechR"))%dopar%{

    ec <- ecs[[i]]
    N2 <- N[i]
    mdl_pars <- unname(unlist(epi_mdl_pars[i,]))

    if (!is.null(priorfunc)) {

      tmp <- optim(mdl_pars, obj_fxn, ec = ec, N = N2, priorfunc = priorfunc,
                   prior = prior[i], timestep = timestep, method = "SANN", control = list(maxit = 1000))

    }

    if (is.null(priorfunc)) {

      tmp <- optim(mdl_pars, obj_fxn, ec = ec, N = N2, timestep = timestep,
                   method = "SANN", control = list(maxit = 1000))

    }

    converged <- tmp$convergence

    curves <- epi_mdl_func(N2, tmp$par[2:length(tmp$par)], timestep, tau + length(ec))

    Kmech <- sum(curves$incidence)

    new_pars <- c(Kmech, tmp$par[2:length(tmp$par)])

    return(list(new_pars, converged))

    gc()

  }

  stopCluster(cl = cluster)

  return(mod_res)
}


#########################
#### Error functions ####
#########################

#' Poisson error term
#'
#' @param ec A numeric vector providing the incidence at each timestep
#' @param pred A numeric vector providing predicted incidence at the each timestep
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#'
#' @returns A numeric object representing the error of the objective function
#' @export
#'
#' @examples
#'
#' poisson_error(ec = observed_cases, pred = model_cases, estK = params[1])
#'
poisson_error <- function(ec, pred, estK) {

  logprob <- sum(dpois(ec, pred, log = TRUE)) +
    dnorm(log10(estK), 0, 1, log = TRUE)

  return(logprob)

}

##

#' Poisson error term penalized towards observed cases
#'
#' @param ec A numeric vector providing the incidence at each timestep
#' @param pred A numeric vector providing predicted incidence at the each timestep. All values must be positive and non-zero.
#' @param estK A numeric object providing the estK parameter used in the mechanistic model. All values must be positive and non-zero.
#'
#' @returns A numeric object representing the error of the objective function
#' @export
#'
#' @examples
#'
#' poisson_error2(ec = observed_cases, pred = model_cases, estK = params[1])
#'
poisson_error2 <- function(ec, pred, estK) {

  logprob <- sum(dpois(ec, pred, log = TRUE)) +
    (dnorm(log10(estK), log10(sum(ec)), 1, log = TRUE)*100)

  return(logprob)

}

##

#' Possion penalty term
#'
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#' @param prior A numeric object providing the prediction from the statistical model
#'
#' @returns A numeric object representing the prior penalty error
#' @export
#'
#' @examples
#'
#' poispen(estK = params[1], prior = stat_result)
#'
poispen <- function(estK, prior) {

  penalty <- dpois(round(estK), prior, log = TRUE)

  return(penalty)

}

##

#' Gaussian error term
#'
#' @param ec A numeric vector providing the incidence at each timestep
#' @param pred A numeric vector providing predicted incidence at the each timestep
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#'
#' @returns A numeric object representing the error of the objective function
#' @export
#'
#' @examples
#'
#' norm_error(ec = observed_cases, pred = model_cases, estK = params[1])
#'
norm_error <- function(ec, pred, estK){

  logprob <- sum(dnorm(ec, pred, 1, log = TRUE)) +
    dnorm(log10(estK), 0, 1, log = TRUE)

  return(logprob)

}

##

#' Gaussian error term penalized towards the observed
#'
#' @param ec A numeric vector providing the incidence at each timestep
#' @param pred A numeric vector providing predicted incidence at the each timestep. All values must be positive and non-zero.
#' @param estK A numeric object providing the estK parameter used in the mechanistic model. All values must be positive and non-zero.
#'
#' @returns A numeric object representing the error of the objective function
#' @export
#'
#' @examples
#'
#' norm_error2(ec = observed_cases, pred = model_cases, estK = params[1])
#'
norm_error2 <- function(ec, pred, estK){

  logprob <- sum(dnorm(ec, pred, 1, log = TRUE)) +
    dnorm(log10(estK), log10(ifelse(sum(ec) > 0, sum(ec), 1)), 1, log = TRUE)

  return(logprob)

}

##

#' Square root penalty term
#'
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#' @param prior A numeric object providing the prediction from the statistical model
#'
#' @returns A numeric object representing the prior penalty error
#' @export
#'
#' @examples
sqrtpen <- function(estK, prior) {
  penalty <- dnorm(sqrt(abs((estK) - prior)), 0, 1, log = TRUE)
  return(penalty)
}

##

#' Diffuse square root penalty term
#'
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#' @param prior A numeric object providing the prediction from the statistical model
#'
#' @returns A numeric object representing the prior penalty error
#' @export
#'
#' @examples
sqrtpen_diffuse <- function(estK, prior) {
  penalty <- dnorm(sqrt(abs((estK) - prior)), 0, 2, log = TRUE)
  return(penalty)
}


##

#' Gaussian penalty term
#'
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#' @param prior A numeric object providing the prediction from the statistical model
#'
#' @returns A numeric object representing the prior penalty error
#' @export
#'
#' @examples
normpen <- function(estK, prior) {
  penalty <- dnorm((abs((estK) - prior)), 0, 1, log = TRUE)
  return(penalty)
}

############################
#### Analysis functions ####
############################

#' Calculate the number of cases averted by allocation of vaccine doses
#'
#' @param method.pred A mathematical expression used to estimate the number of remaining cases
#' @param method.name A character object providing the name of the method
#' @param dat A dataframe providing the location name (`loc`), the true epidemic size (`K`), the observed cases (`observed`), the predicted sizes from the models, and the population size for the location (`N`).
#' @param maxvax A numeric object providing the max number of vaccine doses that can be allocated. The default is `3e6`.
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
#' @import dplyr
#'
#' @export
#'
#' @examples
#'
#' get_infPrevented_beta(method.pred = stat - observed, method.name = "stat",
#'                       dat = my_dataset, maxvax = 3e6)
#'
get_infPrevented_beta <- function(method.pred, method.name, dat = flat_res,
                                   maxvax = 3e6) {
  res <- dat %>%
    ungroup() %>%
    mutate(method = method.name,
           u = runif(n = nrow(dat)),
           est_remaining = {{method.pred}},
           true_remaining = K - observed,
           sortval = ifelse(method == "random", u, (est_remaining)/N)) %>%
    arrange(-sortval) %>%
    mutate(vaxxed = (N) * coverage,
           cumv = cumsum(vaxxed)) %>%
    filter(lag(cumv, n = 1, default = 0) <= maxvax) %>%
    mutate(cumv2 = pmin(cumv, maxvax),
           vaxhere = ifelse(cumv2 == cumv, vaxxed, cumv - cumv2),
           inf_prevented = cumsum((true_remaining*vaxhere/vaxxed) * coverage * VE)) %>%
    select(loc, method, cumv2, inf_prevented)
  return(res)
}

########################
#### Jess functions ####
########################

#' Function for fitting a simple linear statistical model appropriate for
#' passing in to em_func
#'
#' @param x the data frame of covariates to fit based on
#' @param y outcome
#' @return a fit model
#' @export
#'
stat.mdl.lin.fit <- function(x, y) {
  tmp <- data.frame(y = y, x)
  rc <- glm(as.formula(paste0("y ~", paste(names(x), collapse = "+"))), data = tmp, family = "gaussian")
  return(rc)
}

#'Function for doing the predict for doing the predict for the
#' results from stat.mdl.lin.fit
#'
#' @param mdl the fit model
#' @param x covariate matrix
#' @return a vector of predicted final sizes
#' @export
#'
stat.mdl.lin.pred <- function(mdl, x) {
  pred <- pmax(0, predict(mdl, newdata = x))
  return(pred)
}
