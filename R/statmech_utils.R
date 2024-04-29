### Utilities for StatMech Approach ###


##' A function to fit a normal distribution to an epidemic curve
##'
##' Version of the epidemiologic model that fits a normal distribution to an epidemic curve.
##' Initialization fits all parameters, including final size, peak time (i.e., mean) and
##' spread (i.e., standard deviation) (formerly `epimdl.norm.init.beta`)
##'
##' @param ecs the epidemic curves
##'
##' @return a list with final size, peak time and spread
##' @export
##'
epimdl.norm.init <- function(ecs){

  # starting vals
  fs <- pmax(5, 2*sapply(ecs, sum)) #start at twice the observed cases
  peak.time <- sapply(ecs, length) #the consistent choice for peak time is the length of the EC
  spread <- peak.time/3 #the consistent choice for peak time and size
  converged <- rep(0, length(ecs))

  obj.fxn <- function(par,  data){
    fs_i <- par[1]
    peaktime_i <- par[2]
    spread_i <- exp(par[3]) #added exp 9/26/23
    #pred <- fs_i * dnorm(c(1:length(data)), mean = peaktime_i, sd = spread_i)
    pred <- fs_i * (pnorm(1:length(data), peaktime_i, spread_i) -  pnorm((1:length(data))-1, peaktime_i, spread_i))
    logprob <- sum(dpois(data, pred, log = T)) +
      dnorm(log10(fs_i), 0,1, log=TRUE)
    return(-logprob)
  }

  for (i in 1:length(ecs)) {
    ec <- ecs[[i]]
    tmp <- optim(c(fs[i], peak.time[i], log(spread[i])), obj.fxn, data = ec)

    converged[i] <- tmp$convergence
    fs[i] <- (tmp$par[1])
    peak.time[i] <- tmp$par[2]
    spread[i] <- exp(tmp$par[3])

  }

  return(list(K=fs, peak.time = peak.time, spread=spread, converged=converged))
}

##' A function to fit a normal distribution to an epidemic curve with shrinkage towards Kstat
##'
##' Version of the epidemiologic model that fits a normal distribution to an epidemic curve.
##' with shrinkage towards Kstat
##' (previously `epi.mdl.normcurve.fit.gamma`)
##'
##' @param ecs the epidemic curves
##' @param Kstat the values of Kstat
##' @param prev.mdl the previous model, can be null
##' @export
##' @return a list with final size, peak time and spread
##'
epi.mdl.normcurve.fit <- function(ecs, Kstat, prev.mdl=NULL) {

  # starting values
  if (is.null(prev.mdl)) {
    peak.time <- sapply(ecs, length) #the consistent choice for peak time is the length of the EC
    spread <- peak.time/3 #the consistent choice for peak time and size
    fs <- pmax(5, 2*sapply(ecs, sum)) #pmax(1, Kstat) #pmax(1, Kstat) #start at K from stat mdl
  }
  else {
    peak.time <- prev.mdl$peak.time
    spread <- prev.mdl$spread
    fs <- prev.mdl$Kmech
  }

  converged <- rep(0, length(ecs))
  Kmech <- rep(NA, length(ecs))

  ##objective function. POisson error on predicted incidence.
  ob.func <- function(par, ec, priorK) {
    loc.mean <- par[1]
    loc.sd <- exp(par[2])
    loc.K <- exp(par[3]) #force K to be positive
    pred <- loc.K * (pnorm(1:length(ec), loc.mean, loc.sd) -  pnorm((1:length(ec))-1, loc.mean, loc.sd))
    log.prob.epi <- sum(dpois(ec, pred, log=TRUE))
    log.prop.penalty <- dnorm(sqrt(abs(loc.K - priorK)), 0, 1, log=T)
    #log.prob <- log.prob.epi + log.prob.penalty/log.prob.epi
    # can we scale penalty? test out different ones?
    log.prob <- sum(dpois(ec, pred, log=TRUE)) + dpois(round(priorK), loc.K, log = T) #dnorm(log10(abs(loc.K - priorK)), 0, 1, log=T)  #
    # try normal on sqrt scale?

    return(-log.prob)
  }

  for (i in 1:length(ecs)) {
    ec <- ecs[[i]]

    tmp <- optim(c( peak.time[i], log(spread[i]), log(fs[i])), ob.func, ec=ecs[[i]], priorK = Kstat[i])

    converged[i] <- tmp$convergence
    peak.time[i] <- tmp$par[1]
    spread[i] <- exp(tmp$par[2])
    Kmech[i] <- exp(tmp$par[3])

    #print(data.frame(iter, i, Kstat[i], Kmech[i], fs[i], peak.time[i], round(spread[i])))

  }


  return(list(Kmech=Kmech, peak.time = peak.time, spread=spread, converged=converged, ecs=ecs))

}

#
# ##' A function to fit SuperLearner to observed epidemic sizes
# ##'
# ##' Function for doing superlearner fit of the K data appropriate for
# ##' passing in to epiInf.EM
# ##' previously `stat.mdl.sl.fit.beta`
# ##'
# ##' @param x the data frame of covariates to fit based on and outcome K
# ##' @return a fit superlearner model
# ##' @export
#
# stat.mdl.sl.fit <- function(x, y) {
#   require(SuperLearner)
#   require(gam)
#   require(rpart)
#   require(randomForest)
#   rc <- SuperLearner(X = x, Y = y, newX = x, family = "gaussian", SL.library = c("SL.rpart", "SL.randomForest", "SL.glm", 'SL.gam')) #'SL.gam',
#   return(rc)
# }
#
#
# ##' Predict final sizes using SuperLearner model
# ##'
# ##' Function for predicting final epidemic sizes using results from
# ##' stat.mdl.sl.fit
# ##'
# ##' @param mdl the fit model
# ##' @param x dataframe with all variables in `mdl`
# ##'
# ##' @return a vector of predicted final sizes
# ##' @export
# ##'
# stat.mdl.sl.pred <- function(mdl, x) {
#   pred <- pmax(1, predict(mdl, onlySL = T, newx = x)$pred)
#   return(pred)
# }

##' A function to extrapolate from SL model prediction
##'
##' Function for doing the predict for doing the prediction for the
##' results from a stat.mdl.sl.fit (with extrapolation)
##' Use this function for prediction when using the stat model alone.
##'
##' @param mdl the fit model
##' @param x a dataset with the variables in `mdl`
##' @param tau end of follow up (time)
##' @export
##' @return a vector of predicted final sizes

stat.mdl.sl.pred.ex <- function(mdl, x, tau) {
  return(pmax(1, predict(mdl, onlySL = T, newdata = x %>% mutate(t = tau))$pred))
}



##' Use EM to combine results from Stat and Mech models
##'
##' Functoin that fits a general epidemic and statistical model using an EM
##' algorithm.
##'
##' @param epi.curves the individual epidemic curves to fit
##' @param x data frame of covariates.
##' @param starting.K a vector of final size estimates to use as starting values (typically results from an epidemic model alone)
##' @param stat.mdl.fit a function that fits the statistical model with the following parameters and returns a fir model:
##'               - x: a data frame of covariates that must include the final size column K a
##' @param epi.mdl.fit a function that calls the epi model with the following parameter and returns a fit model:
##'               - epi.curves: the epidemic curves.
##'               - K: the projected fornal sizes
##'               - prev.mdl : the previous model. should behave when this is null
##' @param stat.mdl.pred: a function that takes in a statisitcal model and returns a vector of predicted final sizes
##' @param max.iter the maximum number of iterations to run
##' @param threshold the desired precision of the putput
##'
##'
##' @return a vector of estimated final sizes
##' @export

epiInf.EM <- function (epi.curves, x, starting.K,
                          stat.mdl.fit, epi.mdl.fit,
                          stat.mdl.pred,
                          max.iter=100, threshold = 10^-5) {


  iter <- 0 #keep track of iterations

  ##for the first loop through, make sure that we are abover threshold
  iter.diff <- 2*threshold*2

  ##matrices for holding results.
  K <- matrix(nrow=max.iter, ncol= length(epi.curves))
  K.stat <- matrix(nrow=max.iter, ncol= length(epi.curves))
  lastK <- vector(length =nrow(x))

  ##previous epi model
  prev.epi.mdl <- NULL

  ## Loop until ending criteria is met
  while ((iter.diff >= threshold) &
         (iter<max.iter))  {

    iter <- iter+1

    if(iter == 1) lastK <- starting.K
    else lastK <- K[iter-1,]

    ## Fit the statistical model on this iteration
    fit.stat.mdl <-stat.mdl.fit(x, lastK)
    K.stat[iter,] <- stat.mdl.pred(fit.stat.mdl, x)

    ## Fit the epidemic model using the stat stuff
    fit.epi.mdl <- epi.mdl.fit(epi.curves, K.stat[iter,], prev.mdl=prev.epi.mdl)
    prev.epi.mdl <- fit.epi.mdl

    K[iter,] <- fit.epi.mdl$Kmech

    ## update the iter differential if this is the second iteratoin or beyond
    if (iter>1) {
      ##iter.diff <- max(sum(abs(K[iter-1,]-K[iter,])),
      ##                 sum(abs(K.stat[iter-1,]-K.stat[iter,])))
      iter.diff <- sum(abs(K[iter-1,]-K[iter,]))/length(K[iter,])
      #print(cbind(iter, K[iter-1,],K[iter,]))

    }

    ##print(cbind(x$K, K[iter,]))
    ##hist(K.stat[iter,]- K[iter,])
    ##print(which(abs(K.stat[iter,]- K[iter,])>10000))
    #cat(iter, ":",iter.diff ,":", range(K[iter,]),":",
    #    range(K.stat[iter,]),"\n")
  }


  ##return a list with the two Ks and the final Ks
  return(list(K=K[iter,],
              K.stat=K.stat[iter,],
              Khist = K[1:iter,],
              Khist.stat = K.stat[1:iter,],
              epi.parms = fit.epi.mdl))
}


##### STATMECH 2.0 BELOW HERE ######


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
sqrtpen <- function(estK, prior, sd = 1) {
  penalty <- dnorm(sqrt(abs((estK) - prior)), 0, sd, log = TRUE)
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
