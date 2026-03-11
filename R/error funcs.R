## ERROR FXNS ##

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
#'
poisson_error <- function(ec, pred, estK) {

  #pred[which(pred < 1)] <- 1

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

  pred[which(pred < 1)] <- 1

  logprob <- sum(dpois(ec, pred, log = TRUE)) +
    dnorm(log10(estK), log10(sum(ec)), 1, log = TRUE)

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
#'
poispen <- function(estK, prior) {

  prior[which(prior < 1)] <- 1

  penalty <- dpois(round(estK), round(prior), log = TRUE)

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

  if(estK == 0){

    estK <- 1

  }

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

  if(estK == 0){

    estK <- 1

  }

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
#'
#' sqrtpen(estK = pars[1], prior = kstat[1])
#'
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
#'
#' sqrt_diffuse(estK = pars[1], prior = kstat[1])
#'
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
#'
#' normpen(estk = pars[1], prior = kstat[1])
#'
normpen <- function(estK, prior) {

  penalty <- dnorm((abs((estK) - prior)), 0, 1, log = TRUE)

  return(penalty)

}
