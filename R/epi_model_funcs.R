## EPI MODEL FXNS ##

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

#' Optimization of a gaussian epidemic model
#'
#' @param ecs  A list object providing the observed or masked epidemic curves. Each element of the list object is a numeric vector.
#' @param epi_mdl_func A function object providing the function used to simulate the mechanistic model.
#' @param epi_mdl_pars A dataframe providing the parameter values for each parameter (columns) and location (rows).
#' @param error_func A function object providing the error calculation for the optimization of the mechanistic model.
#' @param priorfunc A function object providing the penalty calculation used to penalize the mechanistic model optimization towards the results from the statistical model.
#' @param prior A numeric vector providing the results from the statistical model that the `priorfunc` will penalize towards.
#' @param cores A numeric object providing the number of cores that should be assigned to run the function in parallel. If not running in parallel this can be left as the default `NULL`.
#'
#' @returns A list object consisting of, for each location:
#'            - a dataframe with the updated parameter values
#'            - a numeric object indicating whether the optimization
#'              converged. A `0` indicates convergence.
#'
#' @export
#'
#' @examples
#'
#' fit_norm_model(ecs = my_curves, epi_mdl_func = my_epi_model,
#'                epi_mdl_pars = starting_values, error_func = norm_error,
#'                prior_func = NULL, prior = NULL)
#'
fit_norm_model <- function(ecs, epi_mdl_func, epi_mdl_pars,
                           error_func, priorfunc = NULL, prior = NULL,
                           cores = NULL) {

  obj_fxn <- function(pars, ec, priorfunc = NULL, prior = NULL) {

    pred_curve <- epi_mdl_func(pars, length(ec))

    if(pars[1] == 0){

      pars[1] <- 1

    }

    if(!is.null(priorfunc)) {

      err <- error_func(ec, pred_curve, pars[1]) +
        priorfunc(pars[1], prior)

    }else{

      err <- error_func(ec, pred_curve, pars[1])

    }


    return(-err)

  }


  if(!is.null(cores)){

    cluster <- makeCluster(cores)
    registerDoParallel(cluster)

    mod_res <- foreach(i = seq_along(ecs), .packages = "dplyr")%dopar%{

      ec <- ecs[[i]]

      mdl_pars <- unname(unlist(epi_mdl_pars[i,]))

      if(all(ec == 0)){

        return(list(mdl_pars, NA))

        gc()

      }else{

        if (!is.null(priorfunc)) {

          tmp <- optim(mdl_pars, obj_fxn, ec = ec, priorfunc = priorfunc,
                       prior = prior[i])

        }else{

          tmp <- optim(mdl_pars, obj_fxn, ec = ec)

        }

        converged <- tmp$convergence

        new_pars <- tmp$par

        return(list(new_pars, converged))

        gc()

      }

    }

    stopCluster(cl = cluster)

  }else{

    mod_res <- list()

    for(i in seq_along(ecs)){

      ec <- ecs[[i]]

      mdl_pars <- unname(unlist(epi_mdl_pars[i,]))

      if(all(ec == 0)){

        mod_res[[i]] <- list(mdl_pars, NA)

      }else{

        if (!is.null(priorfunc)) {

          tmp <- optim(mdl_pars, obj_fxn, ec = ec, priorfunc = priorfunc,
                       prior = prior[i])

        }else{

          tmp <- optim(mdl_pars, obj_fxn, ec = ec)

        }

        converged <- tmp$convergence

        new_pars <- tmp$par

        mod_res[[i]] <- list(new_pars, converged)

      }

    }

  }

  return(mod_res)

}

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
#' @param tau A character object providing the number of timesteps that should be simulated for the mechanistic model prediction.
#' @param timestep A numeric object providng the number of days in each time step. For example, a weekly time step would be `timestep = 7`.
#' @param cores A numeric object providing the number of cores that should be assigned to run the function in parallel. If not running in parallel this can be left as the default `NULL`.
#'
#' @returns A list object consisting of, for each location:
#'            - a dataframe with the updated parameter values
#'            - a numeric object indicating whether the optimization
#'              converged. A `0` indicates convergence.
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
                          error_func, priorfunc = NULL, prior = NULL,
                          tau, timestep, cores = NULL) {

  # obj_fxn <- function(pars, ec, N, estK, priorfunc = NULL, prior = NULL, timestep) {

  obj_fxn <- function(pars, ec, N, priorfunc = NULL, prior = NULL, timestep) {

    pred_curve <- epi_mdl_func(N, pars, timestep, (length(ec) - 1))

    estK <- sum(pred_curve$incident)

    if(estK == 0){

      estK <- 1

    }

    if(!is.null(priorfunc)) {

      err <- error_func(ec, pred_curve$incident, estK) +
        priorfunc(estK, prior)

    }else{

      err <- error_func(ec, pred_curve$incident, estK)

    }

    return(-err)

  }

  # obj_fxn <- function(pars, ec, N, priorfunc = NULL, prior = NULL, timestep) {
  #
  #   pred_curve <- epi_mdl_func(N, pars[2:length(pars)], timestep, (length(ec) - 1))
  #
  #   if(pars[1] == 0){
  #
  #     pars[1] <- 1
  #
  #   }
  #
  #   if(!is.null(priorfunc)) {
  #
  #     err <- error_func(ec, pred_curve$incident, pars[1]) +
  #       priorfunc((pars[1]), prior)
  #
  #   }else{
  #
  #     err <- error_func(ec, pred_curve$incident, pars[1])
  #
  #   }
  #
  #   return(-err)
  #
  # }

  if(!is.null(cores)){

    cluster <- makeCluster(cores)
    registerDoParallel(cluster)

    mod_res <- foreach(i = seq_along(ecs), .packages = "dplyr")%dopar%{

      ec <- ecs[[i]]
      N2 <- N[i]
      mdl_pars <- unname(unlist(epi_mdl_pars[i,]))[-1]

      if(all(ec == 0)){

        return(list(c(0, mdl_pars), NA, NA))

        gc()

      }else{

        if (!is.null(priorfunc)) {

          tmp <- optim(mdl_pars, obj_fxn, ec = ec,
                       N = N2, priorfunc = priorfunc,
                       prior = prior[i],
                       timestep = timestep,
                       #method = "SANN",
                       control = list(maxit = 1000))

        }else{

          tmp <- optim(par = mdl_pars, obj_fxn,
                       ec = ec, N = N2,
                       timestep = timestep,
                       #method = "SANN",
                       control = list(maxit = 1000))

        }

        converged <- tmp$convergence

        curves <- epi_mdl_func(N2, tmp$par, timestep, tau + length(ec))$incident

        Kmech <- sum(curves)

        new_pars <- c(Kmech, tmp$par)

        return(list(new_pars, converged, curves))

        gc()

      }

    }

    stopCluster(cl = cluster)

  }else{

    mod_res <- list()

    for(i in seq_along(ecs)){

      ec <- ecs[[i]]
      N2 <- N[i]
      mdl_pars <- unname(unlist(epi_mdl_pars[i,]))[-1]

      if(all(ec == 0)){

        mod_res[[i]] <- list(c(0, mdl_pars), NA, NA)

      }else{

        if (!is.null(priorfunc)) {

          tmp <- optim(mdl_pars, obj_fxn, ec = ec,
                       N = N2, priorfunc = priorfunc,
                       prior = prior[i],
                       timestep = timestep,
                       #method = "SANN",
                       control = list(maxit = 1000))

        }else{

          tmp <- optim(mdl_pars, obj_fxn, ec = ec,
                       N = N2, timestep = timestep,
                       #method = "SANN",
                       control = list(maxit = 1000))

        }

        converged <- tmp$convergence

        curves <- epi_mdl_func(N2, tmp$par, timestep, tau + length(ec))$incident

        Kmech <- sum(curves)

        new_pars <- c(Kmech, tmp$par)

        mod_res[[i]] <- list(new_pars, converged, curves)

      }

    }

  }

  return(mod_res)

}

##
SIR_builder <- function(Sus.to.Inf = list("fitted"),
                        Inf.to.Rec = list("fitted"),
                        detection.prob = list("fixed", 1)){

  function(N, pars, time_step, tau){

    sub_num <- 2

    if(Sus.to.Inf[[1]] == "fixed"){

      beta = Sus.to.Inf[[2]]

    }else{

      beta = pars[[sub_num]]

      sub_num = sub_num + 1

    }

    if(Inf.to.Rec[[1]] == "fixed"){

      phi = Inf.to.Rec[[2]]

    }else{

      phi = pars[[sub_num]]

      sub_num = sub_num + 1

    }

    if(detection.prob[[1]] == "fixed"){

      prob_detect = detection.prob[[2]]

    }else{

      prob_detect = pars[[sub_num]]

      sub_num = sub_num + 1

    }

    I0 = round(N * exp(-abs(pars[1])))

    if(I0 < 1){

      I0 = 1

    }

    cur_state <- c(t = 0, S = N - I0, I = I0 , R = 0,
                   incident = rbinom(1, I0, prob_detect))

    tmp <- list(cur_state)
    i <- 1

    while(i <= round(tau)){

      i = i + 1

      lambda <- (beta * cur_state[3]/N)

      StoI = rbinom(1, cur_state[2], 1-exp(-time_step*lambda))
      StoI = ifelse(is.na(StoI), 0, StoI)

      ItoR = rbinom(1, cur_state[3], 1-exp(-time_step*phi))
      ItoR = ifelse(is.na(ItoR), 0, ItoR)

      cur_state[1] <- cur_state[1] + time_step
      cur_state[2] <- cur_state[2] - StoI
      cur_state[3] <- cur_state[3] + StoI - ItoR
      cur_state[4] <- cur_state[4] + ItoR

      cur_state[5] <- rbinom(1, StoI, prob_detect)

      tmp[[i]] <- cur_state

    }

    epi <- bind_rows(tmp)

    return(epi)

  }

}

####

SEIR_builder <- function(Sus.to.Exp = list("fitted"),
                         Exp.to.Inf = list("fitted"),
                         Inf.to.Rec = list("fitted"),
                         detection.prob = list("fixed", 1)){

  function(N, pars, time_step, tau){

    sub_num <- 2

    if(Sus.to.Exp[[1]] == "fixed"){

      beta = Sus.to.Exp[[2]]

    }else{

      beta = pars[[sub_num]]

      sub_num = sub_num + 1

    }

    if(Exp.to.Inf[[1]] == "fixed"){

      sigma = Exp.to.Inf[[2]]

    }else{

      sigma = pars[[sub_num]]

      sub_num = sub_num + 1

    }

    if(Inf.to.Rec[[1]] == "fixed"){

      phi = Inf.to.Rec[[2]]

    }else{

      phi = pars[[sub_num]]

      sub_num = sub_num + 1

    }

    if(detection.prob[[1]] == "fixed"){

      prob_detect = detection.prob[[2]]

    }else{

      prob_detect = pars[[sub_num]]

      sub_num = sub_num + 1

    }

    I0 = round(N * exp(-abs(pars[1])))

    if(I0 < 1){

      I0 = 1

    }

    cur_state <- c(t = 0, S = N - I0, E = 0, I = I0 , R = 0,
                   incident = rbinom(1, I0, prob_detect))

    tmp <- list(cur_state)
    i <- 1

    while(i <= round(tau)){

      i = i + 1

      lambda <- (beta * cur_state[3]/N)

      StoE = rbinom(1, cur_state[2], 1-exp(-time_step*lambda))
      StoE = ifelse(is.na(StoE), 0, StoE)

      EtoI = rbinom(1, cur_state[3], 1-exp(-time_step*sigma))
      EtoI = ifelse(is.na(EtoI), 0, EtoI)

      ItoR = rbinom(1, cur_state[4], 1-exp(-time_step*phi))
      ItoR = ifelse(is.na(ItoR), 0, ItoR)

      cur_state[1] <- cur_state[1] + time_step
      cur_state[2] <- cur_state[2] - StoE
      cur_state[3] <- cur_state[3] + StoE - EtoI
      cur_state[4] <- cur_state[4] + EtoI - ItoR
      cur_state[5] <- cur_state[5] + ItoR

      cur_state[6] <- rbinom(1, EtoI, prob_detect)

      tmp[[i]] <- cur_state

    }

    epi <- bind_rows(tmp)

    return(epi)

  }

}
