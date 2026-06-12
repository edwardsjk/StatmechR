## EPI MODEL FXNS ##

#' Simulate a an epidemic curve using a Gaussian distribution
#'
#' @description
#' `gaussian_mod()` creates an epidemic curve using a Gaussian distribution, given the following parameters: i) the estimated epidemic size, ii) the estimated peak time of the epidemic, and iii) estimated spread of the epidemic; as well as the length of the epidemic.
#'
#'
#' @param pars A vector providing, in order:
#' - the estimated epidemic size
#' - the estimated peak time of the epidemic
#' - the estimated spread of the epidemic
#' @param times A numeric object providing the number of timesteps that the epidemic should be simulated for
#'
#' @return A vector of incidence values for each timestep
#'
#' @export
#'
#' @examples
#'
#' gaussian_mod(pars = c(250, 15, 5), time = 30)
#'
gaussian_mod <- function(pars,
                         times) {

  if(length(pars) != 3){

    stop("gaussian_mod() requires 3 parameter values")

  }

  fs_i <- unname(unlist(pars[1]))

  peaktime_i <- unname(unlist(pars[2]))

  spread_i <- exp(unname(unlist(pars[3])))

  normpred <- fs_i * (pnorm(1:times, peaktime_i, spread_i) -
                        pnorm((1:times) - 1, peaktime_i, spread_i))

  return(normpred)
}

#' Fit a gaussian curve to epidemic data
#'
#' @description
#' `run_gaussian_model()` fits a gaussian curve to epidemic data to estimate the final epidemic size, using general purpose optimization via `optim()`.
#'
#' @param ecs  A list object providing the observed or masked epidemic curves. Each element of the list object is a numeric vector.
#' @param epi.mdl.pars A dataframe providing the parameter values for each parameter (columns) and location (rows).
#' @param error.func A function object providing the error calculation for the optimization of the mechanistic model.
#' @param prior.func A function object providing the penalty calculation used to penalize the mechanistic model optimization towards the results from the statistical model.
#' @param prior A numeric vector providing the results from the statistical model that the `priorfunc` will penalize towards.
#' @param cores A numeric object providing the number of cores that should be assigned to run the function in parallel. If not running in parallel this can be left as the default `NULL`.
#' @param optim.method A character object providing the optimization method that should be used to fit the epi model to the data. Options are the same as those available in base R's `optim()`.
#' @param optim.control A list object providing additional parameters to control the optimization. Options are the same as those available in base R's `opitm()`.
#'
#' @returns A list object consisting of, for each location:
#'            - a dataframe with the updated parameter values
#'            - a numeric object indicating whether the optimization
#'              converged. A `0` indicates convergence.
#'
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#'
#' @export
#'
#' @examples
#'
#' run_gaussian_model(ecs = my_curves, epi_mdl_pars = starting_values,
#'                    error_func = norm_error, prior_func = NULL,
#'                    prior = NULL)
#'
run_gaussian_model <- function(ecs,
                               epi.mdl.pars,
                               error.func,
                               prior.func = NULL,
                               prior = NULL,
                               cores = NULL,
                               optim.method = "Nelder-Mead",
                               optim.control = NULL) {

  obj_fxn <- function(pars, ec, prior.func = NULL, prior = NULL) {

    pred_curve <- gaussian_mod(pars, length(ec))

    if(!is.null(prior.func)) {

      err <- error.func(ec, pred_curve, pars[1]) +
        prior.func(pars[1], prior)

    }else{

      err <- error.func(ec, pred_curve, pars[1])

    }


    return(-err)

  }


  if(!is.null(cores)){

    cluster <- makeCluster(cores)
    registerDoParallel(cluster)

    mod_res <- foreach(i = seq_along(ecs), .packages = "dplyr")%dopar%{

      ec <- ecs[[i]]

      mdl_pars <- unname(unlist(epi.mdl.pars[i,]))

      if(all(ec == 0)){

        return(list(mdl_pars, NA))

        gc()

      }else{

        if (!is.null(prior.func)) {

          tmp <- optim(mdl_pars, obj_fxn, ec = ec, prior.func = prior.func,
                       prior = prior[i], method = optim.method,
                       control = optim.control)

        }else{

          tmp <- optim(mdl_pars, obj_fxn, ec = ec, method = optim.method,
                       control = optim.control)

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

      mdl_pars <- unname(unlist(epi.mdl.pars[i,]))

      if(all(ec == 0)){

        mod_res[[i]] <- list(mdl_pars, NA)

      }else{

        if (!is.null(prior.func)) {

          tmp <- optim(mdl_pars, obj_fxn, ec = ec, prior.func = prior.func,
                       prior = prior[i], method = optim.method,
                       control = optim.control)

        }else{

          tmp <- optim(mdl_pars, obj_fxn, ec = ec, method = optim.method,
                       control = optim.control)

        }

        converged <- tmp$convergence

        new_pars <- tmp$par

        mod_res[[i]] <- list(new_pars, converged)

      }

    }

  }

  return(mod_res)

}

#' Fit a custom epidemic model to epidemic data
#'
#' @description
#' `run_custom_model()` fits a user defined custom epidemic model to epidemic data to predict final epidemic size, using general purpose optimization via `optim()`.
#'
#'
#' @param ecs  A list object providing the observed or masked epidemic curves. Each element of the list object is a numeric vector.
#' @param N A numeric vector providing the population size in each location.
#' @param epi.mdl.func A function object providing the function used to simulate the mechanistic model.
#' @param epi.mdl.pars A dataframe providing the parameter values for each parameter (columns) and location (rows).
#' @param error.func A function object providing the error calculation for the optimization of the mechanistic model.
#' @param prior.func A function object providing the penalty calculation used to penalize the mechanistic model optimization towards the results from the statistical model.
#' @param prior A numeric vector providing the results from the statistical model that the `prior.func` will penalize towards.
#' @param tau A character object providing the number of timesteps that should be simulated for the mechanistic model prediction.
#' @param timestep A numeric object providng the number of days in each time step. For example, a weekly time step would be `timestep = 7`.
#' @param cores A numeric object providing the number of cores that should be assigned to run the function in parallel. If not running in parallel this can be left as the default `NULL`.
#' @param optim.method A character object providing the optimization method that should be used to fit the epi model to the data. Options are the same as those available in base R's `optim()`.
#' @param optim.control A list object providing additional parameters to control the optimization. Options are the same as those available in base R's `opitm()`.
#'
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
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
#'run_custom_model(ecs = my_curves, epi.mdl.func = my_epi_model,
#'                epi.mdl.pars = starting_values, error.func = norm_error,
#'                prior.func = NULL, prior = NULL, tau = 70,
#'                timestep = 7, cores = NULL, optim.method = "Nelder-Mead",
#'                optim.control = list(maxit = 1000))
#'
run_custom_model <- function(ecs,
                             N,
                             epi.mdl.func,
                             epi.mdl.pars,
                             error.func,
                             prior.func = NULL,
                             prior = NULL,
                             tau,
                             timestep,
                             cores = NULL,
                             optim.method = "Nelder-Mead",
                             optim.control = NULL) {

  obj_fxn <- function(pars, ec, N, prior.func = NULL, prior = NULL, timestep) {

    pred_curve <- epi.mdl.func(N, pars, timestep, (length(ec) + tau))

    estK <- sum(pred_curve$incident[1:length(ec)]) #compare current curve to fitted curve of same length

    estK2 <- sum(pred_curve$incident) #compare prior to total predicted curve

    if (!is.null(prior.func)) {

      err <- error.func(ec, pred_curve$incident, estK) +
        prior.func(estK2, prior)

    } else {

      err <- error.func(ec, pred_curve$incident, estK)

    }

    return(-err)

  }

  if (!is.null(cores)) {

    cluster <- makeCluster(cores)
    registerDoParallel(cluster)

    mod_res <- foreach(i = seq_along(ecs), .packages = "dplyr")%dopar%{

      ec <- ecs[[i]]
      N2 <- N[i]
      mdl_pars <- unname(unlist(epi.mdl.pars[i, ]))

      if (all(ec == 0)) {

        return(list(estK = 0, new_pars = mdl_pars, convergence = NA, curves = NA))

        gc()

      } else {

        if (!is.null(prior.func)) {

          tmp <- optim(mdl_pars, obj_fxn, ec = ec,
                       N = N2, prior.func = prior.func,
                       prior = prior[i],
                       timestep = timestep,
                       method = optim.method,
                       control = optim.control)

        } else {

          tmp <- optim(par = mdl_pars, obj_fxn,
                       ec = ec, N = N2,
                       timestep = timestep,
                       method = optim.method,
                       control = optim.control)

        }

        converged <- tmp$convergence

        curves <- epi.mdl.func(N2, tmp$par, timestep, tau + length(ec))$incident

        Kmech <- sum(curves)

        new_pars <- tmp$par

        return(list(estK = Kmech, params = new_pars,
                    convergence = converged, curves = curves))

        gc()

      }

    }

    stopCluster(cl = cluster)

  }else{

    mod_res <- list()

    for(i in seq_along(ecs)){

      ec <- ecs[[i]]
      N2 <- N[i]
      mdl_pars <- unname(unlist(epi.mdl.pars[i, ]))

      if(all(ec == 0)){

        mod_res[[i]] <- list(estK = 0, params = mdl_pars, convergence = NA, curves = NA)

      }else{

        if (!is.null(prior.func)) {

          tmp <- optim(mdl_pars, obj_fxn, ec = ec,
                       N = N2, prior.func = prior.func,
                       prior = prior[i],
                       timestep = timestep,
                       method = optim.method,
                       control = optim.control)

        }else{

          tmp <- optim(mdl_pars, obj_fxn, ec = ec,
                       N = N2, timestep = timestep,
                       method = optim.method,
                       control = optim.control)

        }

        converged <- tmp$convergence

        curves <- epi.mdl.func(N2, tmp$par, timestep, tau + length(ec))$incident

        Kmech <- sum(curves)

        new_pars <- tmp$par

        mod_res[[i]] <- list(estK = Kmech, params = new_pars,
                             convergence = converged, curves = curves)

      }

    }

  }

  return(mod_res)

}

##
#' Build a simple SIR model
#'
#' @description
#' `build_SIR()` builds a SIR model function that can be passed to `run_custom_model()` and `run_combined_model()`. The SIR model is a simple closed structure that gives the user the ability to fix or fit the parameter values.
#'
#'
#' @param Sus.to.Inf A list object indicating the whether the beta transmission parameter should be `"fitted"` or `"fixed"`. If the beta transmission parameter is `"fixed"`, the second object in the list will be the numeric value of the parameter. If the beta transmission parameter is `"fitted"` no other objects are needed in the list.
#' @param Inf.to.Rec A list object indicating the whether the recovery rate parameter should be `"fitted"` or `"fixed"`. If the recovery rate parameter is `"fixed"`, the second object in the list will be the numeric value of the parameter. If the recovery parameter is `"fitted"` no other objects are needed in the list.
#' @param detection.prob A list object indicating the whether the detection probability parameter should be `"fitted"` or `"fixed"`. If the detection probability parameter is `"fixed"`, the second object in the list will be the numeric value of the parameter. If the detection probability parameter is `"fitted"` no other objects are needed in the list. The default value for the fixed detection probability is 1.
#'
#' @returns A function with the following arguments:
#' -  N = the population size.
#' -  pars = the initial starting values for the fitted parameters.
#' -  time_step = the length of the timestep in the model.
#' -  tau = the length of the simulation
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#'
#' @examples
#'
#' build_SIR(Sus.to.Inf = list("fitted"),
#'           Inf.to.Rec = list("fixed", 0.2),
#'           detection.prob = list("fixed", 0.5))
#'
build_SIR <- function(Sus.to.Inf = list("fitted"),
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

#' Build a simple SEIR model
#'
#' @description
#' `build_SEIR()` builds a SEIR model function that can be passed to `run_custom_model()` and `run_combined_model()`. The SEIR model is a simple closed structure that gives the user the ability to fix or fit the parameter values.
#'
#' @param Sus.to.Exp A list object indicating the whether the beta transmission parameter should be `"fitted"` or `"fixed"`. If the beta transmission parameter is `"fixed"`, the second object in the list will be the numeric value of the parameter. If the beta transmission parameter is `"fitted"` no other objects are needed in the list.
#' @param Exp.to.Inf A list object indicating the whether the progression rate parameter should be `"fitted"` or `"fixed"`. If the progression rate parameter is `"fixed"`, the second object in the list will be the numeric value of the parameter. If the progression rate transmission parameter is `"fitted"` no other objects are needed in the list.
#' @param Inf.to.Rec A list object indicating the whether the recovery rate parameter should be `"fitted"` or `"fixed"`. If the recovery rate parameter is `"fixed"`, the second object in the list will be the numeric value of the parameter. If the recovery rate parameter is `"fitted"` no other objects are needed in the list.
#' @param detection.prob A list object indicating the whether the detection probability parameter should be `"fitted"` or `"fixed"`. If the detection probability parameter is `"fixed"`, the second object in the list will be the numeric value of the parameter. If the detection probability parameter is `"fitted"` no other objects are needed in the list. The default value for the fixed detection probability is 1.
#'
#' @returns A function with the following arguments:
#' -  N = the population size.
#' -  pars = the initial starting values for the fitted parameters.
#' -  time_step = the length of the timestep in the model.
#' -  tau = the length of the simulation
#'
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#'
#' @examples
#'
#' build_SEIR(Sus.to.Exp = list("fitted"),
#'            Exp.to.Inf = list("fixed", 0.33),
#'            Inf.to.Rec = list("fixed", 0.2),
#'            detection.prob = list("fixed", 0.5))
#'
build_SEIR <- function(Sus.to.Exp = list("fitted"),
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
