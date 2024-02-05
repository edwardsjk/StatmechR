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


##' A function to fit SuperLearner to observed epidemic sizes
##'
##' Function for doing superlearner fit of the K data appropriate for
##' passing in to epiInf.EM
##' previously `stat.mdl.sl.fit.beta`
##'
##' @param x the data frame of covariates to fit based on and outcome K
##' @return a fit superlearner model
##' @export

stat.mdl.sl.fit <- function(x, y) {
  require(SuperLearner)
  require(gam)
  require(rpart)
  require(randomForest)
  rc <- SuperLearner(X = x, Y = y, newX = x, family = "gaussian", SL.library = c("SL.rpart", "SL.randomForest", "SL.glm", 'SL.gam')) #'SL.gam',
  return(rc)
}


##' Predict final sizes using SuperLearner model
##'
##' Function for predicting final epidemic sizes using results from
##' stat.mdl.sl.fit
##'
##' @param mdl the fit model
##' @param x dataframe with all variables in `mdl`
##'
##' @return a vector of predicted final sizes
##' @export
##'
stat.mdl.sl.pred <- function(mdl, x) {
  pred <- pmax(1, predict(mdl, onlySL = T, newx = x)$pred)
  return(pred)
}

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
