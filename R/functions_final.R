#### NECESSARY FUNCTIONS ####

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
sim_covariates <- function (N,
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
#' @export
#'
#' @examples
#'
#' #' compute_Rs(covariate_data, b0 = log(2), sd = .1,
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

#' Simulate epidemic data for multiple locations based on covariate data
#'
#' @param N Numeric object providing the number of locations that should be simulated.
#' @param covariates Data frame providing the covariate values for each location.
#' @param pop_max Numeric object providing the maximum population size that could be simulated for a location. The default is `pop_max = 1e+06`
#' @param pop_min Numeric object providing the minimum population size that could be simulated for a location. The default is `pop_min = 1e+03`
#' @param init_I Numeric object providing the number of initial infections per location at the start of the epidemic. The default is `init_I = 10`
#' @param timestep Numeric object providing the time step interval, for example `timestep = 7` would be weekly. The default is `timestep = 7`
#' @param tau Numeric object providing the time frame (in days) that the epidemic should be simulated for.
#' @param offset Numeric object providing the severity of offset applied to simulated outbreaks. Smaller values indicate a larger offset. The default is `offset = 2`.
#'
#' @return A data frame providing the cases of disease for each time step and location, along with the respective covariates that were provided to the function.
#' @export
#'
#'
#' @examples
#
#'  simulate_outbreak(population, covariate_data, pop_max = 1e+06,
#'                    pop_min = 1e+03, init_I = 10, timestep = 7,
#'                    tau = 120, offset = 2)
#'
#'
simulate_outbreak <- function(N, covariates, pop_max = 1000000, pop_min = 1000,
                              init_I = 10, timestep = 7, tau, offset = 2){

  sim_data <- covariates %>%
    mutate(R0 = compute_Rs(covariates), loc = 1 : n(),
           pop = round(runif(n(), pop_min, pop_max)))

  params <- data.frame(I0 = log(sim_data$pop/init_I), beta_a = -log((sim_data$R0*0.1)/2), beta_i = -log(sim_data$R0*0.07))

  epidemics <- data.frame()

  for (i in 1:nrow(sim_data)) {

    epi <- sim_gen_pred(sim_data$pop[i], unlist(params[i,]),
                        time_step = 7, tau)

    epi <- as.data.frame(epi)

    epi$loc <- i
    epi$t <- 1:nrow(epi)

    epidemics <- bind_rows(epidemics, epi)

  }

  offset <- sample(1:round((tau/timestep)/offset), nrow(sim_data), replace = T)

  obs_epis <- epidemics  %>%
    rename(incident = epi) %>%
    #mutate(t = t/timestep)  %>%
    group_by(loc) %>%
    mutate(t = t + offset[loc]) %>%
    #mutate(t = t + sample(1:(tau/timestep)-1, size = n(), replace = T)) %>% # offset dates so they start at different times
    ungroup() %>%
    filter(t <= floor(tau/timestep) & t >= 1) %>%
    select(loc, t, incident)

  obs_full <- data.frame(loc = rep(1:N, each = max(obs_epis$t)),
                         t = rep(1:max(obs_epis$t), times = N),
                         incident = 0)

  obs_epis2 <- left_join(x = rbind(obs_epis,
                                   anti_join(x = obs_full, y = obs_epis,
                                             by = c("t", "loc"))),
                         y = sim_data,
                         by = "loc") %>%
    rename(epiweek_date = t, district = loc, n_cases = incident)

  return(obs_epis2)

}

##

#' Extract the final epidemic size from the epidemic data
#'
#' @param dat A data frame providing the epidemic cases over time for each location
#' @param groups A character object providing the variable for the location level at which the final epidemic size should be calculated
#'
#' @return A data frame providing the location and the corresponding final epidemic size
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

  K <- dat %>%
    ungroup() %>%
    group_by(do.call(pick, grouping)) %>%
    summarise(K = sum(n_cases)) %>%
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
#' @export
#'
#' @examples
#'
#'  plot_true_curves(dat = epidemic_data, X = "epiweek_date",
#'                   plot_group = "district", legend = TRUE)
#'
plot_true_curves <- function(dat, X = "b", plot_group = "a", legend = TRUE){

  grouping <- unlist(lapply(plot_group, function(b){

    which(names(dat) == b)

  }))

  xaxis <- unlist(lapply(X, function(c){

    which(names(dat) == c)

  }))

  plotdat <- dat %>%
    group_by(.[[grouping]], .[[xaxis]]) %>%
    summarize(cases = sum(n_cases)) %>%
    ungroup()

  if(legend == TRUE){

    tplot <- ggplot(plotdat, aes(x = as.Date(`.[[xaxis]]`), y = cases, color = as.factor(`.[[grouping]]`))) +
      geom_line() +
      ylab("New cases") +
      xlab("Week") +
      labs(color = str_to_title(plot_group)) +
      theme(axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18))

  }else{

    tplot <- ggplot(data = plotdat, aes(x = as.Date(`.[[xaxis]]`), y = cases, color = as.factor(`.[[grouping]]`))) +
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
#' @export
#'
#' @examples
#'
#' plot_true_curves(dat = epidemic_data, maskdat = masked_data,
#'                  X = "epiweek_date", plot_group = "district", legend = TRUE)
#'
plot_mask_curves <- function(dat, maskdat, X = "b", plot_group = "a", legend = TRUE){

  grouping1 <- unlist(lapply(plot_group, function(x){

    which(names(maskdat) == x)

  }))

  grouping2 <- unlist(lapply(plot_group, function(x){

    which(names(dat) == x)

  }))

  xaxis <- unlist(lapply(X, function(c){

    which(names(dat) == c)

  }))

  plotdat <- mask_dat %>%
    group_by(.[[grouping1]], .[[xaxis]]) %>%
    summarize(cases = sum(n_cases)) %>%
    ungroup() %>%
    group_by(.[[grouping1]]) %>%
    filter(cases > 0) %>%
    ungroup()

  plotdat_all <- dat %>%
    group_by(.[[grouping2]], .[[xaxis]]) %>%
    summarize(cases = sum(n_cases)) %>%
    ungroup() %>%
    group_by(.[[grouping2]]) %>%
    filter(cases > 0) %>%
    ungroup()

  if(legend == TRUE){

    tplot <- ggplot() +
      geom_line(aes(x = as.Date(`.[[xaxis]]`), y = cases, color = as.factor(`.[[grouping1]]`)),
                data = plotdat) +
      geom_line(aes(x = as.Date(`.[[xaxis]]`), y = cases, color = as.factor(`.[[grouping2]]`)),
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
      geom_line(aes(x = as.Date(`.[[xaxis]]`), y = cases, color = as.factor(`.[[grouping1]]`)),
                data = plotdat) +
      geom_line(aes(x = as.Date(`.[[xaxis]]`), y = cases, color = as.factor(`.[[grouping2]]`)),
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
    mutate(loc = as.numeric(loc)) %>%
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
#' @export
#'
#' @examples
#'
#' stat.mdl.sl.fit.para(x = my_vars, y = epi_size,
#'                      cores = 4, family = "gaussian)
#'
stat.mdl.sl.fit.para <- function(x, y, cores = 1, family = "gaussian") {
  require(SuperLearner)
  require(gam)
  require(rpart)
  require(randomForest)
  require(parallel)

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

##

#' Title
#'
#' @param x
#' @param y
#' @param id
#' @param family
#'
#' @returns
#' @export
#'
#' @examples
stat.mdl.sl.fit <- function(x, y, id, family) {
  require(SuperLearner)
  require(gam)
  require(rpart)
  require(randomForest)

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
#' @export
#'
#' @examples
#'
#' stat.mdl.sl.pred(mdl = SL_model, x = my_vars)
#'
stat.mdl.sl.pred <- function(mdl, x) {

  pred <- pmax(1, predict(mdl, onlySL = T, newdata = x)$pred)

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
#' - `curves` Predicted incidence curves for most recent iteration
#' - `converged` Convergence results from the most recent iteration
#' - `diff` Iteration difference between the last iteration and the one previous
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
                          stat.family, tau, timestep = NULL) {

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

    Kmech2 <- Kmech[iter,]/(population/1000)

    fitstatmdl <- statmdlfit(covdat, Kmech2, cores, stat.family)
    K[iter, ] <- statmdlpred(fitstatmdl, covdat) * (population/1000)

    ## Check iter_diff
    if (iter > 1) {
      iter_diff <- abs(K[iter - 1, ] - K[iter, ])
    }
    #print(rbind(lastK[1:5], Kstat[iter, 1:5], K[iter, 1:5]))
    print(rbind(lastK[1:5], Kmech[iter, 1:5], K[iter, 1:5]))

    gc()

  }

  model_curves <- lapply(fitepimdl, function(x){x[[2]]})
  converged <- lapply(fitepimdl, function(x){x[[3]]})

  return(list(K = K[iter, ],
              #Kstat = Kstat[iter, ],
              Kmech = Kmech[iter, ],
              Khist = K[1:iter, ],
              #Khist.stat = Kstat[1:iter, ],
              Khist.mech = Kmech[1:iter, ],
              epi.parms = fitepimdl,
              param.all = params,
              curves = model_curves,
              converged = converged,
              diff = iter_diff))

}

##

fit_norm_model <- function(ecs, epi_mdl_func, epi_mdl_pars,
                           error_func, priorfunc = NULL, prior = NULL,
                           cores, tau, timestep = NULL) {

#' Title
#'
#' @param pars
#' @param ec
#' @param priorfunc
#' @param prior
#'
#' @returns
#' @export
#'
#' @examples
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

  mod_res <- foreach(i = seq_along(ecs), .packages = "dplyr")%dopar%{

    source("/users/a/b/abagaels/Codes/statmech/R/functions_final.R")

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

    curves <- NA

    return(list(new_pars, curves, converged))

    gc()

  }

  stopCluster(cl = cluster)

  return(mod_res)
}

##

#' Title
#'
#' @param ecs
#' @param N
#' @param epi_mdl_func
#' @param epi_mdl_pred
#' @param epi_mdl_pars
#' @param error_func
#' @param priorfunc
#' @param prior
#' @param cores
#' @param tau
#' @param timestep
#'
#' @returns
#' @export
#'
#' @examples
fit_epi_model <- function(ecs, N, epi_mdl_func, epi_mdl_pred, epi_mdl_pars,
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

  mod_res <- foreach(i = seq_along(ecs), .packages = "dplyr")%dopar%{

    source("/users/a/b/abagaels/Codes/statmech/R/functions_final.R")
    source("/users/a/b/abagaels/Codes/statmech/R/temp_functions.R")

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

    curves <- epi_mdl_pred(N2, tmp$par[2:length(tmp$par)], timestep, tau + length(ec))

    Kmech <- sum(curves)

    new_pars <- c(Kmech, tmp$par[2:length(tmp$par)])

    return(list(new_pars, curves, converged))

    gc()

  }

  stopCluster(cl = cluster)

  return(mod_res)
}


##

#' Title
#'
#' @param ecs
#' @param pop_N
#' @param strt_vals
#' @param errorfxn
#' @param penaltyfunc
#' @param priorval
#' @param cores
#' @param tau
#' @param timestep
#'
#' @returns
#' @export
#'
#' @examples
SEIR_em <- function(ecs, pop_N, strt_vals, errorfxn, penaltyfunc, priorval,
                    cores, tau, timestep) {

  model_fit <- fit_epi_model(ecs, pop_N, SEIR_fit, SEIR_pred, strt_vals,
                             errorfxn, priorfunc = penaltyfunc,
                             prior = priorval, cores, tau, timestep)

  return(model_fit)

}

##

#' Title
#'
#' @param N
#' @param pars
#' @param time_step
#' @param epi_length
#'
#' @returns
#' @export
#'
#' @examples
SEIR_fit <- function(N, pars, time_step, epi_length){

  #pars includes:
  #I0, beta_s, beta_a, beta_w, bac_growth

  # sigma = 0.66
  prop_asymp = 0.8
  # phi_s = 0.2
  # phi_a = 0.2
  prob_detect = 0.425

  I0 = round(N * exp(-abs(pars[1])))

  if(I0 < 1){

    I0 = 1

  }

  beta_s = exp(-abs(pars[2]))
  beta_a = exp(-abs(pars[3]))

  cur_state <- c(t = 0, S = N - I0, E = 0,
                 I_s = round(I0 * (1 - prop_asymp)),
                 I_a = round(I0 * prop_asymp), R = 0,
                 incident = rbinom(1, round(I0 * (1 - prop_asymp)), prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= epi_length){

    i = i + 1

    cur_state <- next_state_SEIR(cur_state, beta_s, beta_a, sigma = 0.66,
                                 prop_asymp, phi_s = 0.2, phi_a = 0.2,
                                 prob_detect, time_step, N)

    tmp[[i]] <- cur_state

  }

  epi <- bind_rows(tmp)

  return(epi)

}

##

#' Title
#'
#' @param N
#' @param pars
#' @param time_step
#' @param tau
#'
#' @returns
#' @export
#'
#' @examples
SEIR_pred <- function(N, pars, time_step, tau){

  #par includes:
  #I0, beta_s, beta_a, beta_w, bac_growth

  # sigma = 0.66
  prop_asymp = 0.8
  # phi_s = 0.2
  # phi_a = 0.2
  prob_detect = 0.425

  I0 = round(N * exp(-abs(pars[1])))

  if(I0 < 1){

    I0 = 1

  }

  beta_s = exp(-abs(pars[2]))
  beta_a = exp(-abs(pars[3]))

  cur_state <- c(t = 0, S = N - I0,
                 E = 0, I_s = round(I0 * (1 - prop_asymp)),
                 I_a = round(I0 * prop_asymp), R = 0,
                 incident = rbinom(1, round(I0 * (1 - prop_asymp)), prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= tau){

    i = i + 1

    cur_state <- next_state_SEIR(cur_state, beta_s, beta_a, sigma = 0.66,
                                 prop_asymp, phi_s = 0.2, phi_a = 0.2,
                                 prob_detect, time_step, N)

    tmp[[i]] <- cur_state

  }

  epi_curve <- bind_rows(tmp)$incident

  return(epi_curve)

}

##

#' Title
#'
#' @param cur_state
#' @param beta_s
#' @param beta_a
#' @param sigma
#' @param prop_asymp
#' @param phi_s
#' @param phi_a
#' @param prob_detect
#' @param time_step
#' @param N
#'
#' @returns
#' @export
#'
#' @examples
next_state_SEIR <- function(cur_state, beta_s, beta_a, sigma,
                            prop_asymp, phi_s, phi_a,
                            prob_detect, time_step,
                            N){

  lambda <- (beta_s * cur_state[4]/N) + (beta_a * cur_state[5]/N)

  StoE = rbinom(1, cur_state[2], 1-exp(-time_step*lambda))
  StoE = ifelse(is.na(StoE), 0, StoE)

  new_infs = rbinom(1, cur_state[3], 1-exp(-time_step*sigma))
  new_infs = ifelse(is.na(new_infs), 0, new_infs)

  IstoR = rbinom(1, cur_state[4], 1-exp(-time_step*phi_s))
  IstoR = ifelse(is.na(IstoR), 0, IstoR)

  IatoR = rbinom(1, cur_state[5], 1-exp(-time_step*phi_a))
  IatoR = ifelse(is.na(IatoR), 0, IatoR)

  EtoIs = round(new_infs * (1-prop_asymp))
  EtoIa = round(new_infs * prop_asymp)

  cur_state[1] <- cur_state[1] + time_step
  cur_state[2] <- cur_state[2] - StoE
  cur_state[3] <- cur_state[3] + StoE - new_infs
  cur_state[4] <- cur_state[4] + EtoIs - IstoR
  cur_state[5] <- cur_state[5] + EtoIa - IatoR
  cur_state[6] <- cur_state[6] + IstoR + IatoR

  cur_state[7] <- rbinom(1, EtoIs, prob_detect)

  return(cur_state)

}

##

#' Title
#'
#' @param ecs
#' @param pop_N
#' @param strt_vals
#' @param errorfxn
#' @param penaltyfunc
#' @param priorval
#' @param cores
#' @param tau
#' @param timestep
#'
#' @returns
#' @export
#'
#' @examples
covid_em <- function(ecs, pop_N, strt_vals, errorfxn, penaltyfunc, priorval,
                     cores, tau, timestep) {

  model_fit <- fit_epi_model(ecs, pop_N, covid_fit, covid_pred, strt_vals,
                             errorfxn, priorfunc = penaltyfunc,
                             prior = priorval, cores, tau, timestep)

  return(model_fit)

}

##

#' Title
#'
#' @param N
#' @param pars
#' @param time_step
#' @param epi_length
#'
#' @returns
#' @export
#'
#' @examples
covid_fit <- function(N, pars, time_step, epi_length){

  #pars includes:
  #I0, beta_p, beta_i, beta_gradient, prob_detect

  # sigma = 0.22
  # phi = 0.36
  #
  # alpha_a = 0.14
  # alpha_i = 0.14

  beta_gradient = exp(-abs(pars[5]))
  prob_detect = exp(-abs(pars[6]))
  prop_asymp = 0.24

  I0 = round(N * exp(-abs(pars[1])))

  if(I0 < 2){

    I0 = 2

  }

  beta_p = exp(-abs(pars[2]))
  beta_i = exp(-abs(pars[3]))
  beta_a = exp(-abs(pars[4]))

  cur_state <- c(t = 0, S = N - I0, E = 0,
                 P = round(I0/2 *(1 - prop_asymp)),
                 A = round(I0/2 * prop_asymp), I = 0, R = 0,
                 incident = rbinom(1, round(I0/2 * (1 - prop_asymp)),
                                   prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= epi_length){

    i = i + 1

    cur_state <- next_state_covid(i, cur_state, beta_p, beta_a, beta_i,
                                  sigma = 0.22, phi = 0.36, alpha_a = 0.14,
                                  alpha_i = 0.14, beta_gradient, prob_detect,
                                  prop_asymp, time_step, N)

    tmp[[i]] <- cur_state

  }

  epi <- bind_rows(tmp)

  return(epi)

}

##

#' Title
#'
#' @param N
#' @param pars
#' @param time_step
#' @param tau
#'
#' @returns
#' @export
#'
#' @examples
covid_pred <- function(N, pars, time_step, tau){

  #pars includes:
  #I0, beta_p, beta_i, beta_gradient, prob_detect

  # sigma = 0.22
  # phi = 0.36
  #
  # alpha_a = 0.14
  # alpha_i = 0.14

  beta_gradient = exp(-abs(pars[5]))
  prob_detect = exp(-abs(pars[6]))
  prop_asymp = 0.24

  I0 = round(N * exp(-abs(pars[1])))

  if(I0 < 1){

    I0 = 2

  }

  beta_p = exp(-abs(pars[2]))
  beta_i = exp(-abs(pars[3]))
  beta_a = exp(-abs(pars[4]))

  cur_state <- c(t = 0, S = N - I0, E = 0,
                 P = round(I0/2 *(1 - prop_asymp)),
                 A = round(I0/2 * prop_asymp), I = 0, R = 0,
                 incident = rbinom(1, round(I0/2 * (1 - prop_asymp)),
                                   prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= tau){

    i = i + 1

    cur_state <- next_state_covid(i, cur_state, beta_p, beta_a, beta_i,
                                  sigma = 0.22, phi = 0.36, alpha_a = 0.14,
                                  alpha_i = 0.14, beta_gradient, prob_detect,
                                  prop_asymp, time_step, N)

    tmp[[i]] <- cur_state

  }

  epi_curve <- bind_rows(tmp)$incident

  return(epi_curve)

}

##

#' Title
#'
#' @param i
#' @param cur_state
#' @param beta_p
#' @param beta_a
#' @param beta_i
#' @param sigma
#' @param phi
#' @param alpha_a
#' @param alpha_i
#' @param beta_gradient
#' @param prob_detect
#' @param prop_asymp
#' @param time_step
#' @param N
#'
#' @returns
#' @export
#'
#' @examples
next_state_covid <- function(i, cur_state, beta_p, beta_a, beta_i, sigma, phi,
                             alpha_a, alpha_i, beta_gradient, prob_detect,
                             prop_asymp, time_step, N){

  if(i >= 24 & i < 93){

    beta_i <- beta_i * 0.19
    beta_p <- beta_p * 0.19
    beta_a <- beta_a * 0.19

  }

  if(i >= 93){

    beta_i_t <- (beta_gradient * i) + (beta_i * 0.19)
    beta_a_t <- (beta_gradient * i) + (beta_a * 0.19)
    beta_p_t <- (beta_gradient * i) + (beta_p * 0.19)

    beta_i <- ifelse(beta_i_t <= beta_i, beta_i_t, beta_i)
    beta_a <- ifelse(beta_a_t <= beta_a, beta_a_t, beta_a)
    beta_p <- ifelse(beta_p_t <= beta_p, beta_p_t, beta_p)

  }

  lambda <- (beta_p * cur_state[4]/N) + (beta_a * cur_state[5]/N) + (beta_i * cur_state[6]/N)

  StoE = rbinom(1, cur_state[2], 1-exp(-time_step*lambda))
  StoE <- ifelse(is.na(StoE) == F, StoE, 0)

  new_infs = rbinom(1, cur_state[3], 1-exp(-time_step*sigma))
  new_infs <- ifelse(is.na(new_infs) == F, new_infs, 0)

  EtoA = round(new_infs * prop_asymp)
  EtoP = round(new_infs * (1 - prop_asymp))

  PtoI = rbinom(1, cur_state[4], 1-exp(-time_step*phi))
  PtoI <- ifelse(is.na(PtoI) == F, PtoI, 0)

  AtoR = rbinom(1, cur_state[5], 1-exp(-time_step*alpha_a))
  AtoR <- ifelse(is.na(AtoR) == F, AtoR, 0)

  ItoR = rbinom(1, cur_state[6], 1-exp(-time_step*alpha_i))
  ItoR <- ifelse(is.na(ItoR) == F, ItoR, 0)

  cur_state[1] <- cur_state[1] + time_step
  cur_state[2] <- cur_state[2] - StoE
  cur_state[3] <- cur_state[3] + StoE - new_infs
  cur_state[4] <- cur_state[4] + EtoP - PtoI
  cur_state[5] <- cur_state[5] + EtoA - AtoR
  cur_state[6] <- cur_state[6] + PtoI - ItoR
  cur_state[7] <- cur_state[7] + AtoR + ItoR
  cur_state[8] <- rbinom(1, PtoI, prob_detect)

  return(cur_state)

}


#########################
#### Error functions ####
#########################

#' Title
#'
#' @param ec
#' @param pred
#' @param estK
#'
#' @returns
#' @export
#'
#' @examples
poisson_error <- function(ec, pred, estK) {

  logprob <- sum(dpois(ec, pred, log = TRUE)) +
    dnorm(log10(estK), 0, 1, log = TRUE)

  return(logprob)

}

##

#' Title
#'
#' @param ec
#' @param pred
#' @param estK
#'
#' @returns
#' @export
#'
#' @examples
poisson_error2 <- function(ec, pred, estK) {

  logprob <- sum(dpois(ec, pred, log = TRUE)) +
    (dnorm(log10(estK), log10(sum(ec)), 1, log = TRUE)*100)

  return(logprob)

}

##

#' Title
#'
#' @param estK
#' @param prior
#'
#' @returns
#' @export
#'
#' @examples
poispen <- function(estK, prior) {

  penalty <- dpois(round(estK), prior, log = TRUE)

  return(penalty)

}

##

#' Title
#'
#' @param ec
#' @param pred
#' @param estK
#'
#' @returns
#' @export
#'
#' @examples
norm_error <- function(ec, pred, estK){

  logprob <- sum(dnorm(ec, pred, 1, log = TRUE)) +
    dnorm(log10(estK), 0, 1, log = TRUE)

  return(logprob)

}

##

#' Title
#'
#' @param ec
#' @param pred
#' @param estK
#'
#' @returns
#' @export
#'
#' @examples
norm_error2 <- function(ec, pred, estK){

  logprob <- sum(dnorm(ec, pred, 1, log = TRUE)) +
    dnorm(log10(estK), log10(ifelse(sum(ec) > 0, sum(ec), 1)), 1, log = TRUE)

  return(logprob)

}

## sqrt penalty
#' Title
#'
#' @param estK
#' @param prior
#'
#' @returns
#' @export
#'
#' @examples
sqrtpen <- function(estK, prior) {
  penalty <- dnorm(sqrt(abs((estK) - prior)), 0, 1, log = TRUE)
  return(penalty)
}


## diffuse sqrt penalty
#' Title
#'
#' @param estK
#' @param prior
#'
#' @returns
#' @export
#'
#' @examples
sqrtpen_diffuse <- function(estK, prior) {
  penalty <- dnorm(sqrt(abs((estK) - prior)), 0, 2, log = TRUE)
  return(penalty)
}


## normal penalty
#' Title
#'
#' @param estK
#' @param prior
#'
#' @returns
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

#' Title
#'
#' @param method.pred
#' @param method.name
#' @param dat
#' @param maxvax
#'
#' @returns
#' @export
#'
#' @examples
get_infPrevented_beta <- function(method.pred, method.name, dat = flat_res,
                                  maxvax = 3e6) {

  #require(tidyverse)
  res <- dat %>%
    ungroup() %>%
    mutate(method = method.name,
           u = runif(n = nrow(dat)),
           est_remaining = {{method.pred}},
           true_remaining = K - observed,
           sortval = ifelse(method == "random", u, (est_remaining)/N)) %>%
    arrange(-sortval) %>%
    #filter(observed>0) %>% # replace with some prob of introduction>0?
    mutate(vaxxed = (N) * coverage,
           cumv = cumsum(vaxxed)) %>%
    filter(cumv<=maxvax) %>%
    mutate(inf_prevented = cumsum((true_remaining) * coverage * VE)) %>%
    select(loc, method, cumv, inf_prevented)

  return(res)

}


#' Title
#'
#' @param method.pred
#' @param method.name
#' @param dat
#' @param maxvax
#'
#' @returns
#' @export
#'
#' @examples
get_infPrevented_beta2 <- function(method.pred, method.name, dat = flat_res,
                                   maxvax = 3e6) {
  #require(tidyverse)
  res <- dat %>%
    ungroup() %>%
    mutate(method = method.name,
           u = runif(n = nrow(dat)),
           est_remaining = {{method.pred}},
           true_remaining = K - observed,
           sortval = ifelse(method == "random", u, (est_remaining)/N)) %>%
    arrange(-sortval) %>%
    # filter(observed>0) %>% # replace with some prob of introduction>0?
    mutate(vaxxed = (N) * coverage,
           cumv = cumsum(vaxxed)) %>%
    filter(lag(cumv, n = 1, default = 0) <= maxvax) %>%
    mutate(cumv2 = pmin(cumv, maxvax),
           vaxhere = ifelse(cumv2 == cumv, vaxxed, cumv - cumv2),
           inf_prevented = cumsum((true_remaining*vaxhere/vaxxed) * coverage * VE)) %>%
    select(loc, method, cumv2, inf_prevented)
  return(res)
}
