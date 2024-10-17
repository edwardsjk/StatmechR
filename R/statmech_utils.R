### Utilities for StatMech Approach ###



##' A function to fit a generic epi model given an epi model and a
##' curve generating function. Currently uses optim with the defualt
##' parameters
##'
##' @param ecs a list of epidemic curves we are trying to fit.
##' @param epi_mdl_func the epidemic model we are trying to fit.
##' Should take a set of parameters and produce a curve of the given length
##' @param epi_mdl_pars a data frame containing the starting parameters
##' for the epidemic model
##' @param error_func the error function we are trying to minimze
##' @param priorfunc prior function used to penalize estimates
##' @param prior a vector of prior values
##' @export
##' @return a data frame with the epidemic model parameters and the covergence result. # nolint

fit_epi_beta <- function(ecs, epi_mdl_func, epi_mdl_pars,
                         error_func, priorfunc = NULL, prior = NULL, ...) {
  ##' takes in the parameter set and the epidemic curve
  obj_fxn <- function(par, ec, priorfunc = NULL, prior = NULL) {
    pred_curve <- epi_mdl_func(par, length(ec))
    if(!is.null(priorfunc)) {
      # err <- sum(error_func(ec, pred_curve, (par[1])),
      #            # stabfunc((par[1])),
      #            priorfunc((par[1]), prior)) #, na.rm = TRUE)
      err <- error_func(ec, pred_curve, (par[1])) +
        priorfunc((par[1]), prior, ...)
    }
    if(is.null(priorfunc)) {
      err <- (error_func(ec, pred_curve, (par[1])))
      #stabfunc((par[1])), na.rm = TRUE)
    }
    return(-err)
  }

  converged <- rep(0, length(ecs))

  for(i in seq_along(ecs)) {
    ec <- ecs[[i]]
    if (!is.null(priorfunc)) {
      tmp <- optim(epi_mdl_pars[i, ], obj_fxn, ec = ec,
                   priorfunc = priorfunc, prior = prior[i])
    }
    if (is.null(priorfunc)) tmp <- optim(epi_mdl_pars[i, ], obj_fxn, ec = ec)

    converged[i] <- tmp$convergence
    epi_mdl_pars[i, ] <- tmp$par
  }

  epi_mdl_pars$converged <- converged

  return(epi_mdl_pars)
}


##' a function to fit a normal model
##' @param par a vector of starting values (final size, peak time, spread)
##' @param times a vector of times to compute predictions
##' @return a vector of predictions at `times`
##' @export
##'
normmdl <- function(par, times) {
  fs_i <- (par[1])
  peaktime_i <- par[2]
  spread_i <- exp(par[3])
  normpred <- fs_i * (pnorm(1:times, peaktime_i, spread_i) -
                        pnorm((1:times) - 1, peaktime_i, spread_i))
  return(normpred)
}

##' A function to compute a poisson error function
##' @param ec a vector of case counts
##' @param pred a vector of predictions
##' @param estK estimated final size
##' @export
##' @return log probability under poisson error structure

poisson_error <- function(ec, pred, estK) {
  logprob <- sum(dpois(ec, pred, log = TRUE)) +
    dnorm(log10(estK), 0, 1, log = TRUE)
  return(logprob)
}

##' a function to stabilize poisson error function
##' built in a and currently NOT USED
##' @export
stabfunc <- function(estK) {
  stab <- dnorm(log10(estK), 0, 1, log = TRUE)
  return(stab)
}

##' a function to make a sqrt penalty
##' @param estK estimated final size
##' @param prior prior estimate of final size
##' @return log probability to be added to error function as a penalty
##' @export
sqrtpen <- function(estK, prior) {
  penalty <- dnorm(sqrt(abs((estK) - prior)), 0, prior/1000, log = TRUE)
  return(penalty)
}


##' a function to make a diffuse sqrt penalty
##' @param estK estimated final size
##' @param prior prior estimate of final size
##' @return log probability to be added to error function as a penalty
##' @export
sqrtpen_diffuse <- function(estK, prior) {
  penalty <- dnorm(sqrt(abs((estK) - prior)), 0, 2, log = TRUE)
  return(penalty)
}


##' a function to make a normal penalty
##' @param estK estimated final size
##' @param prior prior estimate of final size
##' @return log probability to be added to error function as a penalty
##' @export
normpen <- function(estK, prior, sd = 1) {
  penalty <- dnorm((abs((estK) - prior)), 0, sd, log = TRUE)
  return(penalty)
}

##' a function to make a poisson penalty
##' @param estK estimated final size
##' @param prior prior estimate of final size
##' @return log probability to be added to error function as a penalty
##' @export
poispen <- function(estK, prior) {
  penalty <- dpois(round(estK), prior, log = TRUE)
  return(penalty)
}

##' wrapper function to fit a normal epi model within the EM algorithm
##' @param ecs a list of epidemic curves
##' @param strt_vals starting values
##' @param errorfxn error function
##' @param penaltyfunc choice of penalty function
##' @param priorval vector of prior values
##' @return a dataframe with normal model parameters and convergence result
##' @export

normfit_em <- function(ecs, strt_vals, errorfxn, penaltyfunc, priorval) {
  normfit <- fit_epi_beta(ecs, normmdl, strt_vals, poisson_error,
                          priorfunc = penaltyfunc, prior = priorval)
  return(normfit)
}

##' a function to implement EM analysis combining statistical and mechanistic forecasting models
##' @param epi_curves a list of epidemic curves
##' @param covdat a dataframe of covariate data for each geographic area
##' @param initK starting final size values
##' @param epimdlfit name of function used to fit epidemic model
##' @param starting_vals starting values for epidemic model
##' @param error_func error function to use for epi model
##' @param penalty_func penalty function to use for epi model
##' @param statmdlfit function to use to fit stat model
##' @param statmdlpred function to use to get predictions from stat model
##' @param threshold max change between iterations for convergence
##' @param max.iter maximum number of interations
##' @return list of estimated final size K, Khist, Kstat, Kstathist,
##'         and parameters of the epi mdl
##' @export

em_func <- function(epi_curves, covdat, initK,
                    epimdlfit, starting_vals, error_func, penalty_func,
                    statmdlfit, statmdlpred,
                    threshold = 20, max.iter = 100, ...) {

  iter <- 0                   # initialize iter
  iter_diff <- 2 * threshold  # set iter_diff > threshold for first iter

  ##matrices for holding results.
  K <- matrix(nrow = max.iter, ncol = length(epi_curves))
  Kstat <- matrix(nrow = max.iter, ncol = length(epi_curves))
  lastK <- vector(length = nrow(covdat))

  prev_epi_mdl <- starting_vals

  ## Loop until ending criteria is met
  while ((iter_diff >= threshold) & (iter < max.iter)) {
    print(iter)
    iter <- iter + 1                # update iter
    if (iter == 1) lastK <- initK   # set lastK to starting values for first iter
    else lastK <- K[iter - 1, ]     # otherwise, lastK is K from last iter

    ## Fit the statistical model on this iteration
    fitstatmdl <- statmdlfit(covdat, lastK)
    Kstat[iter, ] <- statmdlpred(fitstatmdl, covdat)

    ## Fit the epidemic model using the stat stuff
    fitepimdl <- epimdlfit(epi_curves,
                           strt_vals = data.frame(prev_epi_mdl[, 1], # specific to normal model for now
                                                  prev_epi_mdl[, 2],
                                                  prev_epi_mdl[, 3]),
                           error_func, penalty_func,
                           priorval = Kstat[iter, ], ...)
    prev_epi_mdl <- fitepimdl

    ## Update K
    K[iter, ] <- fitepimdl[, 1]

    ## Check iter_diff
    if (iter > 1) {
      iter_diff <- sum(abs(K[iter - 1, ] - K[iter, ])) / length(K[iter, ])
    }
    print(rbind(lastK[1:5], Kstat[iter, 1:5], K[iter, 1:5]))
  }

  return(list(K = K[iter, ],
              Kstat = Kstat[iter, ],
              Khist = K[1:iter, ],
              Khist.stat = Kstat[1:iter, ],
              epi.parms = fitepimdl))

}


##' Function for doing superlearner fit of the K data appropriate for
##' passing in to em_func
##'
##' @param x the data frame of covariates to fit based on
##' @param y outcome
##' @return a fit superlearner model
##' @export
##'
stat.mdl.sl.fit <- function(x, y) {
  require(SuperLearner)
  require(gam)
  require(rpart)
  require(randomForest)
  rc <- SuperLearner(X = x, Y = y, newX = x, family = "gaussian", SL.library = c("SL.rpart", "SL.randomForest", "SL.glm", 'SL.gam')) #'SL.gam',   
  return(rc)
}


##'Function for doing the predict for doing the predict for the
##' results from a stat.mdl.sl.fit
##'
##' @param mdl the fit model
##'
##' @return a vector of predicted final sizes
##' @export
##'
stat.mdl.sl.pred <- function(mdl, x) {
  pred <- pmax(1, predict(mdl, onlySL = T, newx = x)$pred)
  return(pred)
}

##'Function for doing the predict for doing the predict for the
##' results from a stat.mdl.sl.fit.beta (with extrapolation)
##'
##' @param mdl the fit model
##'
##' @return a vector of predicted final sizes
##' @export
##' 
stat.mdl.sl.pred.ex <- function(mdl, x, tau) {
  return(pmax(1,predict(mdl, onlySL = T, newdata = x %>% mutate(t = tau))$pred))
}


##' Function for fitting a simple linear statistical model appropriate for
##' passing in to em_func
##'
##' @param x the data frame of covariates to fit based on
##' @param y outcome
##' @return a fit model
##' @export
##'
stat.mdl.lin.fit <- function(x, y) {
  tmp <- data.frame(y = y, x)
  rc <- glm(as.formula(paste0("y ~", paste(names(x), collapse = "+"))), data = tmp, family = "gaussian")
  return(rc)
}

##'Function for doing the predict for doing the predict for the
##' results from stat.mdl.lin.fit
##'
##' @param mdl the fit model
##' @param x covariate matrix
##' @return a vector of predicted final sizes
##' @export
##'
stat.mdl.lin.pred <- function(mdl, x) {
  pred <- pmax(1, predict(mdl, newdata = x))
  return(pred)
}
