#### NECESSARY FUNCTIONS ####

###################################
#### Data processing functions ####
###################################

#' Title
#'
#' @param N
#' @param poverty
#' @param climate
#' @param temperature
#' @param water_san
#' @param stand_water
#' @param education
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param covs
#' @param b0
#' @param sd
#' @param b1
#' @param b2
#' @param b3
#' @param temp_scale
#' @param temp_loc
#' @param b4
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param N
#' @param covariates
#' @param pop_max
#' @param pop_min
#' @param init_I
#' @param timestep
#' @param tau
#'
#' @return
#' @export
#'
#' @examples
simulate_outbreak <- function(N, covariates, pop_max = 1000000, pop_min = 1000,
                              init_I = 10, timestep = 7, tau){

  sim_data <- covariates %>%
    mutate(R0 = compute_Rs(covariates), loc = 1 : n(),
           pop = round(runif(n(), pop_min, pop_max)))

  epidemics <- NULL

  for (i in 1:nrow(sim_data)) {

    epi <- fst_slow_OG(sim_data$pop[i], init_I, sim_data$R0[i],
                       time_step = timestep)

    epi$loc <- i
    epidemics <- bind_rows(epidemics, epi)

  }

  offset <- sample(1:round((tau/timestep)/2), nrow(sim_data), replace = T)

  obs_epis <- epidemics  %>%
    mutate(t = t/timestep)  %>%
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

#' Title
#'
#' @param N
#' @param covariates
#' @param pop_max
#' @param pop_min
#' @param init_I
#' @param timestep
#' @param tau
#'
#' @return
#' @export
#'
#' @examples
simulate_outbreak2 <- function(N, covariates, pop_max = 1000000, pop_min = 1000,
                               init_I = 10, timestep = 7, tau){

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

  offset <- sample(1:round((tau/timestep)/2), nrow(sim_data), replace = T)

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

#' Title
#'
#' @param dat
#' @param groups
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param dat
#' @param X
#' @param plot_group
#' @param legend
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param dat
#' @param maskdat
#' @param plot_group
#' @param legend
#'
#' @return
#' @export
#'
#' @examples
plot_mask_curves <- function(dat, maskdat, plot_group = "a", legend = TRUE){

  grouping1 <- unlist(lapply(plot_group, function(x){

    which(names(maskdat) == x)

  }))

  grouping2 <- unlist(lapply(plot_group, function(x){

    which(names(dat) == x)

  }))

  plotdat <- mask_dat %>%
    group_by(.[[grouping1]], epiweek_date) %>%
    summarize(cases = sum(n_cases)) %>%
    ungroup() %>%
    group_by(.[[grouping1]]) %>%
    filter(cases > 0) %>%
    ungroup()

  plotdat_all <- dat %>%
    group_by(.[[grouping2]], epiweek_date) %>%
    summarize(cases = sum(n_cases)) %>%
    ungroup() %>%
    group_by(.[[grouping2]]) %>%
    filter(cases > 0) %>%
    ungroup()

  if(legend == TRUE){

    tplot <- ggplot() +
      geom_line(aes(x = epiweek_date, y = cases, color = as.factor(`.[[grouping1]]`)),
                data = plotdat) +
      geom_line(aes(x = epiweek_date, y = cases, color = as.factor(`.[[grouping2]]`)),
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
      geom_line(aes(x = epiweek_date, y = cases, color = as.factor(`.[[grouping1]]`)),
                data = plotdat) +
      geom_line(aes(x = epiweek_date, y = cases, color = as.factor(`.[[grouping2]]`)),
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

#' Title
#'
#' @param trueK
#' @param iter_results
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param pars
#' @param times
#'
#' @return
#' @export
#'
#' @examples
normmdl <- function(pars, times) {

  fs_i <- unname(unlist(pars[1]))

  peaktime_i <- unname(unlist(pars[2]))

  spread_i <- exp(unname(unlist(pars[3])))

  normpred <- fs_i * (pnorm(1:times, peaktime_i, spread_i) -
                        pnorm((1:times) - 1, peaktime_i, spread_i))

  return(normpred)
}

##

#' Title
#'
#' @param x
#' @param y
#' @param cores
#' @param family
#'
#' @return
#' @export
#'
#' @examples
stat.mdl.sl.fit.para <- function(x, y, cores, family) {
  require(SuperLearner)
  require(gam)
  require(rpart)
  require(randomForest)
  require(parallel)

  cluster = makeCluster(cores)
  clusterEvalQ(cluster, library(SuperLearner))
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

  }

  if(family == "poisson"){

    rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "poisson", cluster = cluster,
                           SL.library = list(c("SL.glm", "screen.glmnet"),
                                             c('SL.gam', "screen.glmnet"),
                                             c("SL.glmnet", "screen.glmnet"),
                                             "SL.glmnet",
                                             c("SL.bartMachine", "screen.glmnet")))

    stopCluster(cluster)

  }

  return(rc)

}

##

#' Title
#'
#' @param mdl
#' @param x
#'
#' @return
#' @export
#'
#' @examples
stat.mdl.sl.pred <- function(mdl, x) {

  pred <- pmax(1, predict(mdl, onlySL = T, newdata = x)$pred)

  return(pred)

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
#' @return
#' @export
#'
#' @examples
norm_em <- function(ecs, pop_N, strt_vals, errorfxn, penaltyfunc, priorval,
                    cores, tau, timestep) {

  model_fit <- fit_norm_model(ecs, normmdl, strt_vals, errorfxn,
                              priorfunc = penaltyfunc, prior = priorval,
                              cores, tau)

  return(model_fit)
}

##

#' Title
#'
#' @param epi_curves
#' @param covdat
#' @param pop_N
#' @param initK
#' @param epimdlfit
#' @param starting_vals
#' @param error_func
#' @param penalty_func
#' @param statmdlfit
#' @param statmdlpred
#' @param threshold
#' @param max.iter
#' @param cores
#' @param stat.family
#' @param tau
#' @param timestep
#'
#' @return
#' @export
#'
#' @examples
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

      lastK <- K[iter - 1, ]} # otherwise, lastK is K from last iter

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

#' Title
#'
#' @param ecs
#' @param epi_mdl_func
#' @param epi_mdl_pars
#' @param error_func
#' @param priorfunc
#' @param prior
#' @param cores
#' @param tau
#' @param timestep
#'
#' @return
#' @export
#'
#' @examples
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
#' @export
#'
#' @examples
covid_fit <- function(N, pars, time_step, epi_length){

  #pars includes:
  #I0, beta_p, beta_i, prob_detect

  # sigma = 0.22
  # phi = 0.36
  #
  # alpha_a = 0.14
  # alpha_i = 0.14

  prob_detect = exp(-abs(pars[4]))
  prop_asymp = 0.24

  I0 = round(N * exp(-abs(pars[1])))

  if(I0 < 2){

    I0 = 2

  }

  beta_p = exp(-abs(pars[2]))
  beta_i = exp(-abs(pars[3]))
  beta_a = beta_i/2

  cur_state <- c(t = 0, S = N - I0, E = 0,
                 P = round(I0/2 *(1 - prop_asymp)),
                 A = round(I0/2 * prop_asymp), I = 0, R = 0,
                 incident = rbinom(1, round(I0/2 * (1 - prop_asymp)),
                                   prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= epi_length){

    i = i + 1

    cur_state <- next_state_covid(cur_state, beta_p, beta_a, beta_i,
                                  sigma = 0.22, phi = 0.36, alpha_a = 0.14,
                                  alpha_i = 0.14, prob_detect, prop_asymp,
                                  time_step, N)

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
#' @return
#' @export
#'
#' @examples
covid_pred <- function(N, pars, time_step, tau){

  #pars includes:
  #I0, beta_p, beta_i, prob_detect

  # sigma = 0.22
  # phi = 0.36
  #
  # alpha_a = 0.14
  # alpha_i = 0.14

  prob_detect = exp(-abs(pars[4]))
  prop_asymp = 0.24

  I0 = round(N * exp(-abs(pars[1])))

  if(I0 < 1){

    I0 = 2

  }

  beta_p = exp(-abs(pars[2]))
  beta_i = exp(-abs(pars[3]))
  beta_a = beta_i/2

  cur_state <- c(t = 0, S = N - I0, E = 0,
                 P = round(I0/2 *(1 - prop_asymp)),
                 A = round(I0/2 * prop_asymp), I = 0, R = 0,
                 incident = rbinom(1, round(I0/2 * (1 - prop_asymp)),
                                   prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= tau){

    i = i + 1

    cur_state <- next_state_covid(cur_state, beta_p, beta_a, beta_i,
                                  sigma = 0.22, phi = 0.36, alpha_a = 0.14,
                                  alpha_i = 0.14, prob_detect,
                                  prop_asymp, time_step, N)

    tmp[[i]] <- cur_state

  }

  epi_curve <- bind_rows(tmp)$incident

  return(epi_curve)

}

##

#' Title
#'
#' @param cur_state
#' @param beta_p
#' @param beta_a
#' @param beta_i
#' @param sigma
#' @param phi
#' @param alpha_a
#' @param alpha_i
#' @param prob_detect
#' @param prop_asymp
#' @param time_step
#' @param N
#'
#' @return
#' @export
#'
#' @examples
next_state_covid <- function(cur_state, beta_p, beta_a, beta_i, sigma, phi,
                             alpha_a, alpha_i, prob_detect,
                             prop_asymp, time_step, N){

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
#' @return
#' @export
#'
#' @examples
sim_gen_em <- function(ecs, pop_N, strt_vals, errorfxn, penaltyfunc, priorval,
                       cores, tau, timestep){

  model_fit <- fit_epi_model(ecs, pop_N, sim_gen_fit, sim_gen_pred, strt_vals,
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
#' @return
#' @export
#'
#' @examples
sim_gen_fit <- function(N, pars, time_step, epi_length){

  #par includes:
  #I0, beta_i, beta_a

  # sigma = 0.1
  prop_asymp = 0.2
  # phi_i = 0.05
  # phi_a = 0.07
  prob_detect = 0.5

  I0 = round(N * exp(-abs(pars[1])))

  if(I0 < 1){

    I0 = 1

  }

  beta_i = exp(-abs(pars[2]))
  beta_a = exp(-abs(pars[3]))

  cur_state <- c(t = 0, S = N - I0,
                 E = 0, I = round(I0 * (1 - prop_asymp)),
                 A = round(I0 * prop_asymp), R = 0,
                 incident = rbinom(1, round(I0 * (1 - prop_asymp)), prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= round(epi_length)){

    i = i + 1

    cur_state <- next_state_gen(cur_state, beta_i, beta_a, sigma = 0.1,
                                prop_asymp, phi_i = 0.07, phi_a = 0.1,
                                prob_detect, time_step, N)

    tmp[[i]] <- cur_state

  }

  epi_curve <- bind_rows(tmp)

  return(epi_curve)

}

##

#' Title
#'
#' @param N
#' @param pars
#' @param time_step
#' @param tau
#'
#' @return
#' @export
#'
#' @examples
sim_gen_pred <- function(N, pars, time_step, tau){

  #par includes:
  #I0, beta_i, beta_a

  # sigma = 0.1
  prop_asymp = 0.2
  # phi_i = 0.05
  # phi_a = 0.07
  prob_detect = 0.5

  I0 = round(N * exp(-abs(pars[1])))

  if(I0 < 1){

    I0 = 1

  }

  beta_a = exp(-abs(pars[2]))
  beta_i = exp(-abs(pars[3]))

  cur_state <- c(t = 0, S = N - I0,
                 E = 0, I = (round(I0 * (1 - prop_asymp))),
                 A = (round(I0 * prop_asymp)), R = 0,
                 incident = rbinom(1, round(I0 * (1 - prop_asymp)), prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= round(tau/time_step)){

    i = i + 1

    cur_state <- next_state_gen(cur_state, beta_i, beta_a, sigma = 0.1,
                                prop_asymp, phi_i = 0.07, phi_a = 0.1,
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
#' @param beta_i
#' @param beta_a
#' @param sigma
#' @param prop_asymp
#' @param phi_i
#' @param phi_a
#' @param prob_detect
#' @param time_step
#' @param N
#'
#' @return
#' @export
#'
#' @examples
next_state_gen <- function(cur_state, beta_i, beta_a, sigma, prop_asymp, phi_i,
                           phi_a, prob_detect, time_step, N){

  lambda <- (beta_i * cur_state[4]/N) + (beta_a * cur_state[5]/N)

  StoE = rbinom(1, cur_state[2], 1-exp(-time_step*lambda))
  StoE = ifelse(is.na(StoE), 0, StoE)

  new_infs = rbinom(1, cur_state[3], 1-exp(-time_step*sigma))
  new_infs = ifelse(is.na(new_infs), 0, new_infs)

  ItoR = rbinom(1, cur_state[4], 1-exp(-time_step*phi_i))
  ItoR = ifelse(is.na(ItoR), 0, ItoR)

  AtoR = rbinom(1, cur_state[5], 1-exp(-time_step*phi_a))
  AtoR = ifelse(is.na(AtoR), 0, AtoR)

  EtoI = round(new_infs * (1-prop_asymp))
  EtoA = round(new_infs * prop_asymp)

  cur_state[1] <- cur_state[1] + time_step
  cur_state[2] <- cur_state[2] - StoE
  cur_state[3] <- cur_state[3] + StoE - new_infs
  cur_state[4] <- cur_state[4] + EtoI - ItoR
  cur_state[5] <- cur_state[5] + EtoA - AtoR
  cur_state[6] <- cur_state[6] + ItoR + AtoR
  cur_state[7] <- rbinom(1, EtoI, prob_detect)

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
#' @return
#' @export
#'
#' @examples
sim_analysis_em <- function(ecs, pop_N, strt_vals, errorfxn, penaltyfunc,
                            priorval, cores, tau, timestep){

  model_fit <- fit_epi_model(ecs, pop_N, sim_analysis_fit, sim_analysis_pred, strt_vals,
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
#' @return
#' @export
#'
#' @examples
sim_analysis_fit <- function(N, pars, time_step, epi_length){

  #par includes:
  #I0, beta_i1, beta_i2

  # sigma = 0.66
  # omega = 0.2
  # phi = 0.2
  prob_detect = 0.425

  I0 = round(N * exp(-abs(pars[1])))

  if(I0 < 1){

    I0 = 1

  }

  beta_i1 = exp(-abs(pars[2]))
  beta_i2 = exp(-abs(pars[3]))

  cur_state <- c(t = 0, S = N - I0,
                 E = 0, I1 = I0,
                 I2 = 0, R = 0,
                 incident = rbinom(1, I0, prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= round(epi_length)){

    i = i + 1

    cur_state <- next_state_analysis(cur_state, beta_i1, beta_i2, sigma = 0.66,
                                     omega = 0.2, phi = 0.2, prob_detect,
                                     time_step, N)

    tmp[[i]] <- cur_state

  }

  epi_curve <- bind_rows(tmp)

  return(epi_curve)

}

##

#' Title
#'
#' @param N
#' @param pars
#' @param time_step
#' @param tau
#'
#' @return
#' @export
#'
#' @examples
sim_analysis_pred <- function(N, pars, time_step, tau){

  #par includes:
  #I0, beta_i1, beta_i2

  # sigma = 0.66
  # omega = 0.2
  # phi = 0.2
  prob_detect = 0.425

  I0 = round(N * exp(-abs(pars[1])))

  if(I0 < 1){

    I0 = 1

  }

  beta_i1 = exp(-abs(pars[2]))
  beta_i2 = exp(-abs(pars[3]))

  cur_state <- c(t = 0, S = N - I0,
                 E = 0, I1 = I0,
                 I2 = 0, R = 0,
                 incident = rbinom(1, I0, prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= round(tau/time_step)){

    i = i + 1

    cur_state <- next_state_analysis(cur_state, beta_i1, beta_i2, sigma = 0.66,
                                     omega = 0.2, phi = 0.2, prob_detect,
                                     time_step, N)

    tmp[[i]] <- cur_state

  }

  epi_curve <- bind_rows(tmp)$incident

  return(epi_curve)

}

##

#' Title
#'
#' @param cur_state
#' @param beta_i1
#' @param beta_i2
#' @param sigma
#' @param omega
#' @param phi
#' @param prob_detect
#' @param time_step
#' @param N
#'
#' @return
#' @export
#'
#' @examples
next_state_analysis <- function(cur_state, beta_i1, beta_i2, sigma, omega, phi,
                                prob_detect, time_step, N){

  lambda <- (beta_i1 * cur_state[4]/N) + (beta_i2 * cur_state[5]/N)

  StoE = rbinom(1, cur_state[2], 1-exp(-time_step*lambda))
  StoE = ifelse(is.na(StoE), 0, StoE)

  EtoI1 = rbinom(1, cur_state[3], 1-exp(-time_step*sigma))
  EtoI1 = ifelse(is.na(EtoI1), 0, EtoI1)

  I1toI2 = rbinom(1, cur_state[4], 1-exp(-time_step*omega))
  I1toI2 = ifelse(is.na(I1toI2), 0, I1toI2)

  I2toR = rbinom(1, cur_state[5], 1-exp(-time_step*phi))
  I2toR = ifelse(is.na(I2toR), 0, I2toR)

  cur_state[1] <- cur_state[1] + time_step
  cur_state[2] <- cur_state[2] - StoE
  cur_state[3] <- cur_state[3] + StoE - EtoI1
  cur_state[4] <- cur_state[4] + EtoI1 - I1toI2
  cur_state[5] <- cur_state[5] + I1toI2 - I2toR
  cur_state[6] <- cur_state[6] + I2toR
  cur_state[7] <- rbinom(1, I1toI2, prob_detect)

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
#' @return
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
#' @return
#' @export
#'
#' @examples
poisson_error2 <- function(ec, pred, estK) {

  if(sum(pred) == 0){

    pred = 0.0000000001

  }

  logprob <- dpois(sum(ec), sum(pred), log = TRUE) +
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
#' @return
#' @export
#'
#' @examples
poisson_error3 <- function(ec, pred, estK) {

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
#' @return
#' @export
#'
#' @examples
poispen <- function(estK, prior) {

  penalty <- dpois(round(estK), prior, log = TRUE)#/mape

  return(penalty)

}

##

#' Title
#'
#' @param ec
#' @param pred
#' @param estK
#'
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @param loc_limit
#'
#' @return
#' @export
#'
#' @examples
get_infPrevented_test <- function(method.pred, method.name, dat = flat_res,
                                  maxvax = 3e6, loc_limit = 50000) {

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
    mutate(vaxxed = ifelse(loc_limit < (N * coverage), loc_limit, N * coverage),
           cumv = cumsum(vaxxed)) %>%
    filter(cumv <= maxvax) %>%
    mutate(inf_prevented = cumsum((true_remaining) * vaxxed/N * VE)) %>%
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
#' @return
#' @export
#'
#' @examples
get_infPrevented2 <- function(method.pred, method.name, dat = flat_res,
                              maxvax = 3e6) {

  #require(tidyverse)
  res <- dat %>%
    ungroup() %>%
    mutate(method = method.name,
           u = runif(n = nrow(dat)),
           est_remaining = {{method.pred}},
           true_remaining = K - observed,
           sortval = ifelse(method == "random", u, est_remaining)) %>%
    arrange(-sortval) %>%
    #filter(observed>0) %>% # replace with some prob of introduction>0?
    mutate(vaxxed = (N) * coverage,
           cumv = cumsum(vaxxed)) %>%
    filter(cumv<=maxvax) %>%
    mutate(inf_prevented = cumsum((true_remaining) * coverage * VE)) %>%
    select(loc, method, cumv, inf_prevented)

  return(res)

}
#######################
#### Old functions ####
#######################
#
# ##
#
# old_em <- function(ecs, pop_N, strt_vals, errorfxn, penaltyfunc, priorval, cores) {
#
#   cholera_fit <- fit_epi_model(ecs, pop_N, SEIR_fit_old, SEIR_pred_old,
#                                strt_vals, errorfxn, priorfunc = penaltyfunc,
#                                prior = priorval, cores)
#
#   return(cholera_fit)
#
# }
#
# ##
#
# SEIR_fit_old <- function(N, par, time_step, epi_length){
#
#   ##par vector/row includes:
#   #I0
#   #R0
#
#   detect_prob = 0.425
#   pct_slow = 0.8
#   sigma = 0.66
#   gamma = 0.2 #length of fast compartment
#   gamma_star = 0.02 #length of slow compartment
#
#   require(tidyverse)
#
#   ## set the transmission parameters - these are daily transmission rates
#   beta <- abs(par[2]) * gamma
#   beta_star <- abs(par[2])/3 * gamma_star/gamma
#
#   ## set the starting state of the epidemic
#   cur_state <- c(t=0, S=N-round(abs(par[1])), E = 0, I=round(abs(par[1]))*0.2,
#                  I_s1=round(abs(par[1]))*0.8, I_s2=0, I_s3=0, R=0,
#                  incident=rbinom(1, round(abs(par[1])), detect_prob))
#
#   tmp <- list(cur_state)
#
#   i <- 1
#
#
#   ## loop through until epidemic length is met
#   while(i <= epi_length){
#     i <- i + 1
#     #update cur state to next state.
#     cur_state <- chol_next_state_old(cur_state, detect_prob, beta, beta_star, pct_slow, sigma, gamma, gamma_star, time_step, N)
#
#     tmp[[i]] <- cur_state
#
#   }
#   # remove bind_rows, but cur_state into list, then bind_rows() the list
#   #add the new current state to the epidemic
#   epi <- bind_rows(tmp)
#
#   return(epi)
# }
#
# ##
#
# SEIR_pred_old <- function(N, par, time_step){
#
#   ##par vector/row includes:
#   #I0
#   #R0
#
#   detect_prob = 0.425
#   pct_slow = 0.8
#   sigma = 0.66
#   gamma = 0.2 #length of fast compartment
#   gamma_star = 0.02 #length of slow compartment
#
#   require(tidyverse)
#
#   ## set the transmission parameters - these are daily transmission rates
#   beta <- abs(par[2]) * gamma
#   beta_star <- abs(par[2])/3 * gamma/gamma_star
#
#   ## set the starting state of the epidemic
#   cur_state <- c(t=0, S=N-round(abs(par[1])), E=0, I=round(abs(par[1]))*0.2,
#                  I_s1=round(abs(par[1]))*0.8, I_s2=0, I_s3=0, R=0,
#                  incident=rbinom(1, round(abs(par[1])), detect_prob))
#   tmp <- list(cur_state)
#   i <- 1
#
#
#   while(i < 504/7){
#
#     i <- i + 1
#
#     #update cur state to next state.
#     cur_state <- chol_next_state_old(cur_state, detect_prob, beta, beta_star, pct_slow, sigma, gamma, gamma_star, time_step, N)
#
#     tmp[[i]] <- cur_state
#
#   }
#
#   epi_curve <- bind_rows(tmp)[,9]
#
#   return(epi_curve)
#
# }
#
# ##
#
# next_state_old <- function(cur_state,
#                            detect_prob,
#                            beta,
#                            beta_star,
#                            pct_slow,
#                            sigma,
#                            gamma,
#                            gamma_star,
#                            time_step,
#                            N) {
#
#   require(tidyverse)
#
#   ##calculate the FOI
#   lambda <- beta * cur_state[4]/N + beta_star * (cur_state[5] +
#                                                    cur_state[6] + cur_state[7])/N
#
#   ##calculate all the transitions
#
#   StoE = rbinom(1, cur_state[2], 1-exp(-time_step*lambda))
#   nw_infs = rbinom(1, cur_state[3], 1-exp(-time_step*sigma))
#   EtoI_s1 = rbinom(1, nw_infs, pct_slow)
#   ItoR = rbinom(1, cur_state[4], 1-exp(-time_step*gamma))
#   I_s1toI_s2 = rbinom(1, cur_state[5], 1-exp(-time_step*gamma_star))
#   I_s2toI_s3 = rbinom(1, cur_state[6], 1-exp(-time_step*gamma_star))
#   I_s3toR = rbinom(1, cur_state[7], 1-exp(-time_step*gamma_star))
#
#   ls <- list(StoE, nw_infs, EtoI_s1, ItoR, I_s1toI_s2, I_s2toI_s3, I_s3toR)
#
#   transitions <- unlist(lapply(ls, function(x){
#
#     if(is.na(x)){
#
#       x = 0
#
#     }else{
#
#       x = x
#
#     }
#
#   }))
#
#   EtoI = transitions[2] -  transitions[3]
#
#   cur_state[1] <- cur_state[1] + time_step
#   cur_state[2] <-  cur_state[2] - transitions[1]
#   cur_state[3] <- cur_state[3]+ transitions[1] - EtoI - transitions[3]
#   cur_state[4] <-  cur_state[4] + EtoI - transitions[4]
#   cur_state[5] <-  cur_state[5] + transitions[3] - transitions[5]
#   cur_state[6] <- cur_state[6] + transitions[5] - transitions[6]
#   cur_state[7] <-  cur_state[7] + transitions[6] - transitions[7]
#   cur_state[8] <-  cur_state[8] + transitions[4] + transitions[7]
#   cur_state[9] <-  rbinom(1, transitions[2], detect_prob)
#
#   return(cur_state)
# }
#
# ##
#
# stat.mdl.sl.fit.CV <- function(x, y, family) {
#   require(SuperLearner)
#   require(gam)
#   require(rpart)
#   require(randomForest)
#   require(parallel)
#
#   # cluster = makeCluster(cores)
#   # clusterEvalQ(cluster, library(SuperLearner))
#   # clusterSetRNGStream(cluster, 1)
#
#   if(family == "gaussian"){
#
#     rc <- CV.SuperLearner(X = x, Y = y, V = 10, family = "gaussian",
#                           #cluster = cluster,
#                           parallel = "multicore",
#                           SL.library = list(c("SL.rpart", "screen.randomForest"),
#                                             c("SL.randomForest", "screen.randomForest"),
#                                             c("SL.glm", "screen.randomForest"),
#                                             c('SL.gam', "screen.randomForest"),
#                                             c("SL.glmnet", "screen.randomForest"),
#                                             c("SL.xgboost", "screen.randomForest"),
#                                             c("SL.svm", "screen.randomForest"),
#                                             "SL.glmnet",
#                                             c("SL.bartMachine", "screen.randomForest")))
#
#     #stopCluster(cluster)
#
#   }
#
#   if(family == "poisson"){
#
#     rc <- CV.SuperLearner(X = x, Y = y, V = 10, family = "poisson",
#                           #cluster = cluster,
#                           parallel = "multicore",
#                           SL.library = list(c("SL.glm", "screen.randomForest"),
#                                             c('SL.gam', "screen.randomForest"),
#                                             c("SL.glmnet", "screen.randomForest"),
#                                             "SL.glmnet",
#                                             c("SL.bartMachine", "screen.randomForest")))
#
#     #stopCluster(cluster)
#
#   }
#
#   return(rc)
#
# }
#
# ##
#
# SEIRB_em <- function(ecs, pop_N, strt_vals, errorfxn, penaltyfunc, priorval, cores, tau, timestep, mape = NULL) {
#
#   cholera_fit <- fit_epi_model(ecs, pop_N, SEIRB_fit, SEIRB_pred, strt_vals,
#                                errorfxn, priorfunc = penaltyfunc,
#                                prior = priorval, cores, tau, timestep)
#
#   return(cholera_fit)
#
# }
#
# ##
#
# SEIRB_fit <- function(N, pars, time_step, epi_length){
#
#   #pars includes:
#   #I0, beta_s, beta_a, beta_w, bac_growth
#
#   sigma = 0.66
#   prop_asymp = 0.8
#   phi_s = 0.2
#   phi_a = 0.2
#   bac_decay = 0.033
#   omega_s = 10
#   omega_a = 1
#   envir_0.5 = 10^6
#   prob_detect = 0.425
#   bac_growth = abs(pars[5])
#
#   I0 = round(N * exp(-abs(pars[1])))
#
#   if(I0 < 1){
#
#     I0 = 1
#
#   }
#
#   beta_s = exp(-abs(pars[2]))
#   beta_a = exp(-abs(pars[3]))
#   beta_w = exp(-abs(pars[4]))
#
#   cur_state <- c(t = 0, S = N - I0, E = 0,
#                  I_s = round(I0 * (1 - prop_asymp)),
#                  I_a = round(I0 * prop_asymp), R = 0, W = 1,
#                  incident = rbinom(1, round(I0 * (1 - prop_asymp)), prob_detect))
#
#   tmp <- list(cur_state)
#   i <- 1
#
#   while(i <= epi_length){
#
#     i = i + 1
#
#     cur_state <- next_state_SEIRB(cur_state, beta_s, beta_a, beta_w, sigma, prop_asymp, phi_s, phi_a, bac_growth, bac_decay, omega_s, omega_a, envir_0.5, prob_detect, time_step, N)
#
#     tmp[[i]] <- cur_state
#
#   }
#
#   epi <- bind_rows(tmp)
#
#   return(epi)
#
# }
#
# ##
#
# SEIRB_pred <- function(N, pars, time_step, tau){
#
#   #par includes:
#   #I0, beta_s, beta_a, beta_w, bac_growth
#
#   sigma = 0.66
#   prop_asymp = 0.8
#   phi_s = 0.2
#   phi_a = 0.2
#   bac_decay = 0.033
#   omega_s = 10
#   omega_a = 1
#   envir_0.5 = 10^6
#   prob_detect = 0.425
#   bac_growth = abs(pars[5])
#
#   I0 = round(N * exp(-abs(pars[1])))
#
#   if(I0 < 1){
#
#     I0 = 1
#
#   }
#
#   beta_s = exp(-abs(pars[2]))
#   beta_a = exp(-abs(pars[3]))
#   beta_w = exp(-abs(pars[4]))
#
#   cur_state <- c(t = 0, S = N - I0,
#                  E = 0, I_s = round(I0 * (1 - prop_asymp)),
#                  I_a = round(I0 * prop_asymp), R = 0, W = 1,
#                  incident = rbinom(1, round(I0 * (1 - prop_asymp)), prob_detect))
#
#   tmp <- list(cur_state)
#   i <- 1
#
#   while(i <= tau){
#
#     i = i + 1
#
#     cur_state <- next_state_SEIRB(cur_state, beta_s, beta_a, beta_w, sigma, prop_asymp, phi_s, phi_a, bac_growth, bac_decay, omega_s, omega_a, envir_0.5, prob_detect, time_step, N)
#
#     tmp[[i]] <- cur_state
#
#   }
#
#   epi_curve <- bind_rows(tmp)$incident
#
#   return(epi_curve)
#
# }
#
# ##
#
# next_state_SEIRB <- function(cur_state, beta_s, beta_a, beta_w, sigma,
#                              prop_asymp, phi_s, phi_a, bac_growth, bac_decay,
#                              omega_s, omega_a, envir_0.5, prob_detect, time_step,
#                              N){
#
#   lambda <- (beta_s * cur_state[4]/N) + (beta_a * cur_state[5]/N) + ((beta_w * cur_state[7])/(envir_0.5 + cur_state[7]))
#
#   StoE = rbinom(1, cur_state[2], 1-exp(-time_step*lambda))
#   new_infs = rbinom(1, cur_state[3], 1-exp(-time_step*sigma))
#   IstoR = rbinom(1, cur_state[4], 1-exp(-time_step*phi_s))
#   IatoR = rbinom(1, cur_state[5], 1-exp(-time_step*phi_a))
#   Wdecay = rbinom(1, cur_state[7], 1-exp(-time_step*bac_decay))
#
#   ls <- list(StoE, new_infs, IstoR, IatoR, Wdecay)
#
#   transitions <- unlist(lapply(ls, function(x){
#
#     if(is.na(x)){
#
#       x = 0
#
#     }else{
#
#       x = x
#
#     }
#
#   }))
#
#   EtoIs = round(transitions[2] * (1-prop_asymp))
#   EtoIa = round(transitions[2] * prop_asymp)
#   shedIs = cur_state[4] * (omega_s * time_step)
#   shedIa = cur_state[5] * (omega_a * time_step)
#
#   cur_state[1] <- cur_state[1] + time_step
#   cur_state[2] <- cur_state[2] - transitions[1]
#   cur_state[3] <- cur_state[3] + transitions[1] - transitions[2]
#   cur_state[4] <- cur_state[4] + EtoIs - transitions[3]
#   cur_state[5] <- cur_state[5] + EtoIa - transitions[4]
#   cur_state[6] <- cur_state[6] + transitions[3] + transitions[4]
#
#   W_change <- ((bac_growth*7*cur_state[7]) - transitions[5]) + shedIs + shedIa
#
#   if(W_change >= 0){
#
#     cur_state[7] <- W_change
#
#   }else{
#
#     cur_state[7] <- 0
#
#   }
#
#   cur_state[8] <- rbinom(1, cur_state[4], prob_detect)
#
#   return(cur_state)
#
# }
#
# fst_slow_OG <- function (N, I0, R0,
#                          detect_prob = 0.33,
#                          pct_slow = 0.5,
#                          sigma = 1/1.5,
#                          gamma = 1/2,
#                          gamma_star = 1/20,
#                          time_step = 0.25) {
#
#   require(tidyverse)
#
#   ## set the transmission parameters
#   beta <- R0 * gamma
#   beta_star <- R0/3 * gamma_star/gamma
#
#   cur_state <- c(t=0, S=N-I0, E=0, I=I0, I_s1=0, I_s2=0, I_s3=0, R=0,
#                  incident=rbinom(1,I0,detect_prob))
#   tmp <- list(cur_state)
#   i <- 1
#
#   while((cur_state["E"] + cur_state["I"] + cur_state["I_s1"] + cur_state["I_s2"] + cur_state["I_s3"]) > 0) {
#     i <- i + 1
#     #update cur state to next state.
#     cur_state <- next_state_OG(cur_state, detect_prob, beta, beta_star, pct_slow, sigma, gamma, gamma_star, time_step, N)
#
#     tmp[[i]] <- cur_state
#
#   }
#   # remove bind_rows, but cur_state into list, then bind_rows() the list
#   #add the new current state to the epidemic
#   epi <- bind_rows(tmp)
#
#   return(epi)
# }
#
# ##
#
# next_state_OG <- function( cur_state,
#                            detect_prob,
#                            beta,
#                            beta_star,
#                            pct_slow,
#                            sigma,
#                            gamma,
#                            gamma_star,
#                            time_step,
#                            N) {
#
#   require(tidyverse)
#
#   ##calculate the FOI
#   lambda <- beta * cur_state["I"]/N + beta_star * (cur_state["I_s1"] +
#                                                      cur_state["I_s2"] + cur_state["I_s3"])/N
#
#   ##calculate all the transitions
#   StoE <- rbinom(1, cur_state["S"], 1-exp(-time_step*lambda))
#   nw_infs <- rbinom(1, cur_state["E"], 1-exp(-time_step*sigma))
#   EtoI_s1<- rbinom(1, nw_infs, pct_slow)
#   EtoI <- nw_infs -  EtoI_s1
#   ItoR <- rbinom(1, cur_state["I"], 1-exp(-time_step*gamma))
#   I_s1toI_s2 <- rbinom(1, cur_state["I_s1"], 1-exp(-time_step*gamma_star))
#   I_s2toI_s3 <- rbinom(1, cur_state["I_s2"], 1-exp(-time_step*gamma_star))
#   I_s3toR <- rbinom(1, cur_state["I_s3"], 1-exp(-time_step*gamma_star))
#
#   cur_state["t"] <- cur_state["t"]+time_step
#   cur_state["S"] <-  cur_state["S"]-StoE
#   cur_state["E"] <- cur_state["E"]+ StoE - EtoI - EtoI_s1
#   cur_state["I"] <-  cur_state["I"] + EtoI - ItoR
#   cur_state["I_s1"] <-  cur_state["I_s1"] + EtoI_s1 - I_s1toI_s2
#   cur_state["I_s2"] <- cur_state["I_s2"] + I_s1toI_s2 - I_s2toI_s3
#   cur_state["I_s3"] <-  cur_state["I_s3"] + I_s2toI_s3 - I_s3toR
#   cur_state["R"] <-  cur_state["R"] + ItoR + I_s3toR
#   cur_state["incident"] <-  rbinom(1,nw_infs,detect_prob)
#
#   return(cur_state)
# }
#
# ##
