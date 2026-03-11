## COMBINED FXNS ##

#' Combined mechanistic and statistic iterative model
#'
#' @param epi_curves A list object providing the observed or masked epidemic curves. Each element of the list object is a numeric vector.
#' @param covdat A data frame object providing the covariate values for each location.
#' @param pop_N A numeric vector providing the population size for each location.
#' @param stat.model A character string indicating the statistical methodology that should be used. The options are a SuperLearner `"SL"` (the default) or a linear regression `"lin"`.
#' @param epi.model A character string indicating the epidemic model that should be used.The options are a Gaussian model `"gaussian"` or a custom model `"custom"`.
#' @param custom_model A function object providing the custom epidemic model if applicable. If using the Gaussian model, this argument can be left as the default `NULL`.
#' @param starting_vals A data frame providing the starting parameter values for the mechanistic model.
#' @param error_func A function object providing the error calculation for the optimization of the mechanistic model.
#' @param penalty_func A function object providing the penalty calculation used to penalize the mechanistic model optimization towards the results from the statistical model.
#' @param tau A character object providing the number of timesteps that should be simulated for the mechanistic model prediction. If using the gaussian model, allow tau to be the default `NULL`.
#' @param timestep A numeric object providing the number of days in each time step. For example, a weekly time step would be `timestep = 7`.If using the gaussian model, allow timestep to be the default `NULL`.
#' @param stat.family A character object providing the error distribution family for the statistical model. The options are `"gaussian"` and `"poisson"`.
#' @param threshold A numeric object providing a threshold for the iteration difference. The iteration difference is difference between the current iteration prediction and previous iteration prediction. Once this drops below the threshold among all locations, the iterative process will stop. The default is `20`
#' @param max.iter A numeric object providing the maximum number of iterations that should be completed. If the threshold has not been met by this number of iterations, the iterative process will stop. The default is `100`
#' @param epi.parallel A logical object indicating if the epidemic model should be run in parallel. The default is `FALSE`.
#' @param stat.parallel A logical object indicating if the statistical model should be run in parallel. The default is `FALSE`
#' @param cores A numeric object providing the number of cores that should be assigned to run the function in parallel. If not running in parallel this can be left as the default `NULL`.
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
#'
#'
#' @export
#'
#' @examples
#'
#' em_func_model(epi_curves = my_curves, covdat = my_covariates,
#'               stat.model = "SL", epi_model = "gaussian",
#'               custom_model = NULL, starting_vals = my_params,
#'               error_func = poisson_error2, penaltyfunc = poispen,
#'               tau = NULL, timestep = NULL, stat.family = "gaussian",
#'               threshold = 5, max.iter = 100, epi.parallel = TRUE,
#'               stat.parallel = FALSE, cores = 4)
#'
#'
#'
em_func_model <- function(epi_curves, covdat, pop_N,
                          stat.model = "SL", epi.model = "gaussian",
                          library.SL = NULL,
                          custom_model = NULL, starting_vals, error_func,
                          penaltyfunc, tau = NULL, timestep = NULL,
                          stat.family = "gaussian", threshold_pct = 0.05,
                          max.iter = 100, epi.parallel = F,
                          stat.parallel = F, cores = NULL){

  obs <- unname(unlist(lapply(epi_curves, function(x){sum(x)})))

  ##matrices for holding results.
  K <- matrix(nrow = max.iter, ncol = length(epi_curves))
  Kmech <- matrix(nrow = max.iter, ncol = length(epi_curves))

  prev_epi_mdl <- starting_vals

  threshold <- pmax(median(obs * threshold_pct), 1)

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

          fitepimdl <- fit_norm_model(ecs = epi_curves,
                                      epi_mdl_func = normmdl,
                                      epi_mdl_pars = starting_vals,
                                      error_func = error_func,
                                      priorfunc = NULL,
                                      prior = NULL)

        }else{

          fitepimdl <- fit_norm_model(ecs = epi_curves,
                                      epi_mdl_func = normmdl,
                                      epi_mdl_pars = starting_vals,
                                      error_func = error_func,
                                      priorfunc = NULL,
                                      prior = NULL,
                                      cores = cores)

        }

      }

      if(epi.model == "custom"){

        #prev_epi_mdl[, 2:ncol(prev_epi_mdl)] <- log(prev_epi_mdl[, -1])

        if(epi.parallel == FALSE){

          fitepimdl <- fit_epi_model(ecs = epi_curves,
                                     N = pop_N,
                                     epi_mdl_func = custom_model,
                                     epi_mdl_pars = starting_vals,
                                     error_func = error_func,
                                     priorfunc = NULL,
                                     prior = NULL,
                                     tau = tau,
                                     timestep = timestep)

        }else{

          fitepimdl <- fit_epi_model(ecs = epi_curves,
                                     N = pop_N,
                                     epi_mdl_func = custom_model,
                                     epi_mdl_pars = starting_vals,
                                     error_func = error_func,
                                     priorfunc = NULL,
                                     prior = NULL,
                                     cores = cores,
                                     tau = tau,
                                     timestep = timestep)

        }

      }

    }else{

      if(epi.model == "gaussian"){

        if(epi.parallel == FALSE){

          fitepimdl <- fit_norm_model(ecs = epi_curves,
                                      epi_mdl_func = normmdl,
                                      epi_mdl_pars = starting_vals,
                                      error_func = error_func,
                                      priorfunc = penaltyfunc,
                                      prior = lastK)

        }else{

          fitepimdl <- fit_norm_model(ecs = epi_curves,
                                      epi_mdl_func = normmdl,
                                      epi_mdl_pars = starting_vals,
                                      error_func = error_func,
                                      priorfunc = penaltyfunc,
                                      prior = lastK,
                                      cores = cores)

        }

      }

      if(epi.model == "custom"){

        #prev_epi_mdl[, 2:ncol(prev_epi_mdl)] <- log(prev_epi_mdl[, -1])

        if(epi.parallel == FALSE){

          fitepimdl <- fit_epi_model(ecs = epi_curves,
                                     N = pop_N,
                                     epi_mdl_func = custom_model,
                                     epi_mdl_pars = starting_vals,
                                     error_func = error_func,
                                     priorfunc = penaltyfunc,
                                     prior = lastK,
                                     tau = tau,
                                     timestep = timestep)

        }else{

          fitepimdl <- fit_epi_model(ecs = epi_curves,
                                     N = pop_N,
                                     epi_mdl_func = custom_model,
                                     epi_mdl_pars = starting_vals,
                                     error_func = error_func,
                                     priorfunc = penaltyfunc,
                                     prior = lastK,
                                     cores = cores,
                                     tau = tau,
                                     timestep = timestep)

        }

      }

    }

    tmp <- lapply(fitepimdl, function(x){x[[1]]})

    prev_epi_mdl <- as.data.frame(do.call(rbind, tmp))

    # if(epi.model == "custom"){
    #
    #   prev_epi_mdl[, 2:ncol(prev_epi_mdl)] <- exp(abs(prev_epi_mdl[, -1]))
    #
    # }

    Kmech[iter, ] <- prev_epi_mdl[, 1]

    ## Fit the statistical model on this iteration

    Kmech2 <- Kmech[iter,]/(pop_N/1000)

    if(stat.model == "SL"){

      if(stat.parallel == FALSE){

        fitstatmdl <- stat.mdl.sl.fit(x = covdat, y = Kmech2,
                                      family = stat.family,
                                      library.SL = library.SL)

        K[iter, ] <- stat.mdl.sl.pred(fitstatmdl, covdat) * (pop_N/1000)


      }else{

        fitstatmdl <- stat.mdl.sl.fit.para(x = covdat, y = Kmech2,
                                           cores = cores, family = stat.family,
                                           library.SL = library.SL)

        K[iter, ] <- stat.mdl.sl.pred(fitstatmdl, covdat) * (pop_N/1000)

      }

    }

    if(stat.model == "linear"){

      fitstatmdl <- stat.mdl.lin.fit(x = covdat, y = Kmech2)

      K[iter, ] <- stat.mdl.lin.pred(fitstatmdl, covdat) * (pop_N/1000)

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
