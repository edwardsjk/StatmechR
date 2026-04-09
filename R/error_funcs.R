## ERROR FXNS ##

#' Poisson error function penalizing towards 0
#'
#' @description
#' `poisson_error()` provides an error values based on two terms: i) the sum of the distance between the observed incidence of the epidemic and the predicted incidence at each timestep; and ii) the distance between the total predicted incidence and 0. This error is then used by the optimization method to fit the model to the observed data.
#'
#' @param ec A numeric vector providing the incidence at each timestep
#' @param pred A numeric vector providing predicted incidence at the each timestep
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#'
#' @returns A numeric object representing the error of the objective function
#'
#' @export
#'
#' @examples
#'
#' poisson_error(ec = observed_cases, pred = model_cases,
#'               estK = sum(model_cases))
#'
#'
poisson_error <- function(ec,
                          pred,
                          estK) {

  rlang::check_required(ec)
  rlang::check_required(pred)
  rlang::check_required(estK)

  check_NAs(ec, thresold = "any")
  check_NAs(pred, threshold = "any")
  check_NAs(estK, threshold = "any")

  check_empty(ec, "vector"); check_empty(pred, "vector")
  check_empty(estK, "vector")

  check_class(ec, "numeric"); check_class(pred, "numeric")
  check_class(estK, "numeric")

  check_positive0(estK); check_length(estK, 1)

  logprob <- sum(dpois(ec, pred, log = TRUE)) +
    dnorm(log10(estK), 0, 1, log = TRUE)

  return(logprob)

}

##

#' Poisson error function penalizing towards observed cases
#'
#' #' @description
#' `poisson_error2()` provides an error values based on two terms: i) the sum of the distance between the observed incidence of the epidemic and the predicted incidence at each timestep; and ii) the distance between the total predicted incidence and the total observed cases. This error is then used by the optimization method to fit the model to the observed data.
#'
#' @param ec A numeric vector providing the incidence at each timestep
#' @param pred A numeric vector providing predicted incidence at the each timestep. All values must be positive and non-zero.
#' @param estK A numeric object providing the estK parameter used in the mechanistic model. All values must be positive and non-zero.
#'
#' @returns A numeric object representing the error of the objective function
#'
#' @export
#'
#' @examples
#'
#' poisson_error2(ec = observed_cases, pred = model_cases, estK = sum(model_cases))
#'
poisson_error2 <- function(ec,
                           pred,
                           estK) {

  rlang::check_required(ec)
  rlang::check_required(pred)
  rlang::check_required(estK)

  check_NAs(ec, thresold = "any")
  check_NAs(pred, threshold = "any")
  check_NAs(estK, threshold = "any")

  check_empty(ec, "vector"); check_empty(pred, "vector")
  check_empty(estK, "vector")

  check_class(ec, "numeric"); check_class(pred, "numeric")
  check_class(estK, "numeric")

  check_positive0(estK); check_length(estK, 1)

  logprob <- sum(dpois(ec, pred, log = TRUE)) +
    dnorm(log10(estK), log10(sum(ec)), 1, log = TRUE)

  return(logprob)

}

##

#' Poisson penalty function penalizing towards the statistical prediction
#'
#' @description
#' `poispen()` provides a penalty for the optimization based on the distance between the total predicted incidence from the epi model and the predicted incidence for the stat model.This penalty value is combined with the value produced by the error functions during the optimization.
#'
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#' @param prior A numeric object providing the prediction from the statistical model
#'
#' @returns A numeric object representing the prior penalty error
#' @export
#'
#' @examples
#'
#' poispen(estK = sum(model_cases), prior = stat_cases)
#'
#'
poispen <- function(estK,
                    prior) {

  rlang::check_required(estK)
  rlang::check_required(prior)

  check_NAs(prior, threshold = "any")
  check_NAs(estK, threshold = "any")
  check_empty(prior, "vector"); check_empty(estK, "vector")
  check_class(prior, "numeric"); check_class(estK, "numeric")
  check_positive0(prior); check_positive0(estK)
  check_length(prior, 1); check_length(estK, 1)

  penalty <- dpois(round(estK), round(prior), log = TRUE)

  return(penalty)

}

##

#' Gaussian error function penalizing towards 0
#'
#' @description
#' `gaussian_error()` provides an error values based on two terms: i) the sum of the distance between the observed incidence of the epidemic and the predicted incidence at each timestep; and ii) the distance between the total predicted incidence and 0. This error is then used by the optimization method to fit the model to the observed data.
#'
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
#' gaussian_error(ec = observed_cases, pred = model_cases,
#'                estK = sum(model_cases))
#'
gaussian_error <- function(ec,
                           pred,
                           estK){

  rlang::check_required(ec)
  rlang::check_required(pred)
  rlang::check_required(estK)

  check_NAs(ec, thresold = "any")
  check_NAs(pred, threshold = "any")
  check_NAs(estK, threshold = "any")

  check_empty(ec, "vector"); check_empty(pred, "vector")
  check_empty(estK, "vector")

  check_class(ec, "numeric"); check_class(pred, "numeric")
  check_class(estK, "numeric")

  check_positive0(estK); check_length(estK, 1)

  logprob <- sum(dnorm(ec, pred, 1, log = TRUE)) +
    dnorm(log10(estK), 0, 1, log = TRUE)

  return(logprob)

}

##

#' Gaussian error function penalizing towards the observed cases
#'
#' @description
#' `gaussian_error2()` provides an error values based on two terms: i) the sum of the distance between the observed incidence of the epidemic and the predicted incidence at each timestep; and ii) the distance between the total predicted incidence and total observed cases. This error is then used by the optimization method to fit the model to the observed data.
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
#' gaussian_error2(ec = observed_cases, pred = model_cases, estK = sum(model_cases))
#'
gaussian_error2 <- function(ec,
                            pred,
                            estK){

  rlang::check_required(ec)
  rlang::check_required(pred)
  rlang::check_required(estK)

  check_NAs(ec, thresold = "any")
  check_NAs(pred, threshold = "any")
  check_NAs(estK, threshold = "any")

  check_empty(ec, "vector"); check_empty(pred, "vector")
  check_empty(estK, "vector")

  check_class(ec, "numeric"); check_class(pred, "numeric")
  check_class(estK, "numeric")

  check_positive0(estK); check_length(estK, 1)

  logprob <- sum(dnorm(ec, pred, 1, log = TRUE)) +
    dnorm(log10(estK), log10(ifelse(sum(ec) > 0, sum(ec), 1)), 1, log = TRUE)

  return(logprob)

}

##

#' Square root penalty function penalizing towards the statistical prediction
#'
#' @description
#' `sqrtpen()` provides a penalty for the optimization based on the distance between the total predicted incidence from the epi model and the predicted incidence for the stat model.This penalty value is combined with the value produced by the error functions during the optimization.
#'
#'
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#' @param prior A numeric object providing the prediction from the statistical model
#'
#' @returns A numeric object representing the prior penalty error
#'
#' @export
#'
#' @examples
#'
#' sqrtpen(estK = sum(model_cases), prior = stat_cases)
#'
sqrtpen <- function(estK,
                    prior) {

  rlang::check_required(estK)
  rlang::check_required(prior)

  check_NAs(prior, threshold = "any")
  check_NAs(estK, threshold = "any")
  check_empty(prior, "vector"); check_empty(estK, "vector")
  check_class(prior, "numeric"); check_class(estK, "numeric")
  check_positive0(prior); check_positive0(estK)
  check_length(prior, 1); check_length(estK, 1)

  penalty <- dnorm(sqrt(abs((estK) - prior)), 0, 1, log = TRUE)

  return(penalty)

}

##

#' Square root penalty function penalizing towards the statistical prediction
#'
#' @description
#' `sqrtpen_diff()` provides a penalty for the optimization based on the distance between the total predicted incidence from the epi model and the predicted incidence for the stat model.This penalty value is combined with the value produced by the error functions during the optimization.
#'
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#' @param prior A numeric object providing the prediction from the statistical model
#'
#' @returns A numeric object representing the prior penalty error
#'
#' @export
#'
#' @examples
#'
#' sqrtpen_diff(estK = sum(model_cases), prior = stat_cases)
#'
sqrtpen_diff <- function(estK,
                         prior) {

  rlang::check_required(estK)
  rlang::check_required(prior)

  check_NAs(prior, threshold = "any")
  check_NAs(estK, threshold = "any")
  check_empty(prior, "vector"); check_empty(estK, "vector")
  check_class(prior, "numeric"); check_class(estK, "numeric")
  check_positive0(prior); check_positive0(estK)
  check_length(prior, 1); check_length(estK, 1)

  penalty <- dnorm(sqrt(abs((estK) - prior)), 0, 2, log = TRUE)

  return(penalty)

}

##

#' Gaussian penalty function penalizing towards the statistical prediction
#'
#' @description
#' `gausspen()` provides a penalty for the optimization based on the distance between the total predicted incidence from the epi model and the predicted incidence for the stat model.This penalty value is combined with the value produced by the error functions during the optimization.
#'
#' @param estK A numeric object providing the estK parameter used in the mechanistic model prediction
#' @param prior A numeric object providing the prediction from the statistical model
#'
#' @returns A numeric object representing the prior penalty error
#' @export
#'
#' @examples
#'
#' gausspen(estk = sum(model_cases), prior = stat_cases)
#'
gausspen <- function(estK,
                     prior) {

  rlang::check_required(estK)
  rlang::check_required(prior)

  check_NAs(prior, threshold = "any")
  check_NAs(estK, threshold = "any")
  check_empty(prior, "vector"); check_empty(estK, "vector")
  check_class(prior, "numeric"); check_class(estK, "numeric")
  check_positive0(prior); check_positive0(estK)
  check_length(prior, 1); check_length(estK, 1)

  penalty <- dnorm((abs((estK) - prior)), 0, 1, log = TRUE)

  return(penalty)

}
