## COMBINED FXNS ##

#' Run the combined mechanistic-statistical EpiStat model
#'
#' @description
#' `run_combined_model()` implements the EpiStat model framework - a framework combining mechanistic and statistical models to improve prediction of final epidemic size. The outputs from each model type are used to inform the next in an iterative fashion until a maximum nubmer of iterations is met or the outputs converge.
#'
#' @param epi.curves A list object providing the observed or masked epidemic curves. Each element of the list object is a numeric vector.
#' @param covdat A data frame object providing the covariate values for each location.
#' @param pop.N A numeric vector providing the population size for each location.
#' @param stat.model A character string providing the statistical methodology that should be used. The options are a SuperLearner `"SL"` (the default) or a linear regression `"linear"`.
#' @param library.SL A list object providing the SL algorithms and screeners to be used in the Superlearner. Using the option `NULL` which will use the default set of algorithms and screeners.
#' @param epi.model A character string providing the epidemic model that should be used.The options are a Gaussian model `"gaussian"` or a custom model `"custom"`.
#' @param custom.model A function object providing the custom epidemic model if applicable. If using the Gaussian model, this argument can be left as the default `NULL`.
#' @param starting.vals A data frame providing the starting parameter values for the mechanistic model.
#' @param optim.method A character object providing the optimization method that should be used to fit the epi model to the data. Options are the same as those available in base R's `optim()`.
#' @param optim.control A list object providing additional parameters to control the optimization. Options are the same as those available in base R's `opitm()`.
#' @param error.func A function object providing the error calculation for the optimization of the mechanistic model.
#' @param penalty.func A function object providing the penalty calculation used to penalize the mechanistic model optimization towards the results from the statistical model.
#' @param tau A character object providing the number of timesteps that should be simulated for the mechanistic model prediction. If using the gaussian model, allow tau to be the default `NULL`.
#' @param timestep A numeric object providing the number of days in each time step. For example, a weekly time step would be `timestep = 7`.If using the gaussian model, allow timestep to be the default `NULL`.
#' @param stat.family A character object providing the error distribution family for the statistical model. The options are `"gaussian"` and `"poisson"`.
#' @param threshold.pct A numeric object providing the threshold that must be met for the model to stop. This given as a proportion (i.e. between 0 and 1) and relates to the proportion of the observed cases for each location. The value is compared to the difference between the epi and stat model components at each iteration. Once the difference falls below the proportion the model/framework will stop.
#' @param max.iter A numeric object providing the maximum number of iterations that should be completed. If the threshold has not been met by this number of iterations, the iterative process will stop. The default is `100`
#' @param epi.parallel A logical object indicating if the epidemic model should be run in parallel. The default is `FALSE`.
#' @param stat.parallel A logical object indicating if the statistical model should be run in parallel. The default is `FALSE`
#' @param cores A numeric object providing the number of cores that should be assigned to run the function in parallel. If not running in parallel this can be left as the default `NULL`.
#'
#'
#'
#' @return A list object containing the following objects:
#' - `K` Most recent prediction of epidemic size from combined model
#' - `Kmech` Most recent prediction of epidemic size from the mechanistic component
#' - `Khist` Epidemic size predictions from the combined model at each iteration
#' - `Khist.mech` Epidemic size predictions from mechanistic component at each iteration
#' - `epi.params` Most recently optimized model parameters
#' - `params` Optimized model parameters for the last epidemic iteration
#' - `converged` Convergence results from the most recent iteration
#' - `diff` Iteration difference between the last iteration and the one previous
#'
#' @export
#'
#' @examples
#'
#' run_combined_model(epi.curves = my_curves, covdat = my_data,
#'                    pop.N = population, stat.model = "SL",
#'                    library.SL = NULL, epi.model = "gaussian",
#'                    custom.model = NULL, starting.vals = my_vals,
#'                    optim.method = "Nelder-Mead",
#'                    optim.control = list(maxit = 1000),
#'                    error.func = poisson_error, penalty.func = poispen,
#'                    tau = 60, timestep = 7, stat.family = "gaussian",
#'                    threshold.pct = 0.05, max.iter = 100, epi.parallel = F,
#'                    stat.parallel = F, cores = NULL)
#'
#'
run_combined_model <- function(epi.curves,
                               covdat,
                               pop.N,
                               stat.model = "SL",
                               library.SL = NULL,
                               epi.model = "gaussian",
                               custom.model = NULL,
                               starting.vals,
                               optim.method = "Nelder-Mead",
                               optim.control = NULL,
                               error.func,
                               penalty.func,
                               tau = NULL,
                               timestep = NULL,
                               stat.family = "gaussian",
                               threshold.pct = 0.05,
                               max.iter = 100,
                               epi.parallel = F,
                               stat.parallel = F,
                               cores = NULL){

  obs <- unname(unlist(lapply(epi.curves, function(x){sum(x)})))

  ##matrices for holding results.
  K <- matrix(nrow = max.iter, ncol = length(epi.curves))
  Kmech <- matrix(nrow = max.iter, ncol = length(epi.curves))

  prev_epi_mdl <- starting.vals

  threshold <- pmax(median(obs * threshold.pct), 1)

  iter <- 0                   # initialize iter
  iter_diff <- 2 * threshold  # set iter_diff > threshold for first iter

  ## Loop until ending criteria is met
  while (any(iter_diff >= threshold) & (iter < max.iter)) {
    print(iter)
    iter <- iter + 1                # update iter
    if(iter == 1){

      lastK <- prev_epi_mdl[,1]

    }else{  # set lastK to starting values for first iter

      lastK <- K[iter - 1, ]

    } # otherwise, lastK is K from last iter

    ## Force lastK and first parameter to be >= obs

    lastK[which(lastK < obs)] <- obs[which(lastK < obs)]

    prev_epi_mdl[which(prev_epi_mdl[,1] < obs), 1] <- obs[which(prev_epi_mdl[,1] < obs)]

    if(iter == 1){

      if(epi.model == "gaussian"){

        if(epi.parallel == FALSE){

          fitepimdl <- run_gaussian_model(ecs = epi.curves,
                                          epi.mdl.func = gaussian_mod,
                                          epi.mdl.pars = starting.vals,
                                          error.func = error.func,
                                          prior.func = NULL,
                                          prior = NULL,
                                          optim.method,
                                          optim.control)

        }else{

          fitepimdl <- run_gaussian_model(ecs = epi.curves,
                                          epi.mdl.func = gaussian_mod,
                                          epi.mdl.pars = starting.vals,
                                          error.func = error.func,
                                          prior.func = NULL,
                                          prior = NULL,
                                          cores = cores,
                                          optim.method,
                                          optim.control)

        }

      }

      if(epi.model == "custom"){

        #prev_epi_mdl[, 2:ncol(prev_epi_mdl)] <- log(prev_epi_mdl[, -1])

        if(epi.parallel == FALSE){

          fitepimdl <- run_custom_model(ecs = epi.curves,
                                        N = pop.N,
                                        epi.mdl.func = custom.model,
                                        epi.mdl.pars = starting.vals,
                                        error.func = error.func,
                                        prior.func = NULL,
                                        prior = NULL,
                                        tau = tau,
                                        timestep = timestep,
                                        optim.method,
                                        optim.control)

        }else{

          fitepimdl <- run_custom_model(ecs = epi.curves,
                                        N = pop.N,
                                        epi.mdl.func = custom.model,
                                        epi.mdl.pars = starting.vals,
                                        error.func = error.func,
                                        prior.func = NULL,
                                        prior = NULL,
                                        cores = cores,
                                        tau = tau,
                                        timestep = timestep,
                                        optim.method,
                                        optim.control)

        }

      }

    }else{

      if(epi.model == "gaussian"){

        if(epi.parallel == FALSE){

          fitepimdl <- run_gaussian_model(ecs = epi.curves,
                                          epi.mdl.func = gaussian_mod,
                                          epi.mdl.pars = starting.vals,
                                          error.func = error.func,
                                          prior.func = penalty.func,
                                          prior = lastK,
                                          optim.method,
                                          optim.control)

        }else{

          fitepimdl <- run_gaussian_model(ecs = epi.curves,
                                          epi.mdl.func = gaussian_mod,
                                          epi.mdl.pars = starting.vals,
                                          error.func = error.func,
                                          prior.func = penalty.func,
                                          prior = lastK,
                                          cores = cores,
                                          optim.method,
                                          optim.control)

        }

      }

      if(epi.model == "custom"){

        #prev_epi_mdl[, 2:ncol(prev_epi_mdl)] <- log(prev_epi_mdl[, -1])

        if(epi.parallel == FALSE){

          fitepimdl <- run_custom_model(ecs = epi.curves,
                                        N = pop.N,
                                        epi.mdl.func = custom.model,
                                        epi.mdl.pars = starting.vals,
                                        error.func = error.func,
                                        prior.func = penalty.func,
                                        prior = lastK,
                                        tau = tau,
                                        timestep = timestep,
                                        optim.method,
                                        optim.control)

        }else{

          fitepimdl <- run_custom_model(ecs = epi.curves,
                                        N = pop.N,
                                        epi.mdl.func = custom.model,
                                        epi.mdl.pars = starting.vals,
                                        error.func = error.func,
                                        prior.func = penalty.func,
                                        prior = lastK,
                                        cores = cores,
                                        tau = tau,
                                        timestep = timestep,
                                        optim.method,
                                        optim.control)

        }

      }

    }

    tmp <- lapply(fitepimdl, function(x) {x[[1]]})

    prev_epi_mdl <- as.data.frame(do.call(rbind, tmp))

    # if(epi.model == "custom"){
    #
    #   prev_epi_mdl[, 2:ncol(prev_epi_mdl)] <- exp(abs(prev_epi_mdl[, -1]))
    #
    # }

    Kmech[iter, ] <- prev_epi_mdl[, 1]

    ## Fit the statistical model on this iteration

    Kmech2 <- Kmech[iter,] / (pop_N / 1000)

    if(stat.model == "SL"){

      if(stat.parallel == FALSE){

        fitstatmdl <- SL_fit(x = covdat, y = Kmech2,
                             family = stat.family,
                             library.SL = library.SL)

        K[iter, ] <- SL_pred(fitstatmdl, covdat) * (pop_N / 1000)


      }else{

        fitstatmdl <- SL_fit(x = covdat, y = Kmech2,
                             cores = cores, family = stat.family,
                             library.SL = library.SL)

        K[iter, ] <- SL_pred(fitstatmdl, covdat) * (pop_N / 1000)

      }

    }

    if(stat.model == "linear"){

      fitstatmdl <- linear_fit(x = covdat, y = Kmech2)

      K[iter, ] <- linear_pred(fitstatmdl, covdat) * (pop_N / 1000)

    }

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
              param.all = prev_epi_mdl,
              converged = converged,
              diff = iter_diff))

}
