## STAT MODEL FXNS ##

#' Train a SuperLearner ensemble
#'
#' @description
#' `SL_fit()` trains a SuperLearner machine learning ensemble on outbreak data. There are options to define the library of algorithms and functionality to run in parallel using multiple cores.
#'
#' @param x A dataframe providing the variables for each location.
#' @param y A numeric vector providing the outcome for each location.
#' @param family A character object providing the error distribution. Currenly only `gaussian` (the default) and `binomial` are supported.
#' @param library.SL A list object providing the SL algorithms and screeners to be used in the Superlearner. Using the option `NULL` which will use the default set of algorithms and screeners.
#' @param parallel A logical object indicating whether the SuperLearner should be run in parallel. The default is `FALSE`.
#' @param cores A numeric object providing the number of cores to be used in the parallel process. The default value is `NULL`.
#'
#' @return A `SuperLearner` object containing the trained model
#'
#' @import SuperLearner
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterSetRNGStream
#' @importFrom parallel stopCluster
#'
#' @export
#'
#' @examples
#'
#' SL_fit(x = my_vars, y = epi_size, family = "gaussian",
#'            library.SL = my_library, parallel = F, cores = NULL)
#'
#'

SL_fit <- function(x,
                    y,
                    family = "gaussian",
                    CV.control = list(V = 10),
                    library.SL = NULL,
                    parallel = F,
                    cores = NULL) {

  if(parallel == T){

    cluster = makeCluster(cores)
    clusterEvalQ(cluster,
                 {library(SuperLearner)
                   options(mc.cores = 1)})
    clusterSetRNGStream(cluster, 1)

    if(is.null(library.SL)){

      if(family == "gaussian"){

        rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "gaussian",
                               cluster = cluster, cvControl = CV.control,
                               SL.library = list(c("SL.rpart", "All"),
                                                 c("SL.glm", "All"),
                                                 c("SL.gam", "All"),
                                                 c("SL.svm", "All"),
                                                 c("SL.glmnet", "All")))

        stopCluster(cluster)
        gc()

      }

      if(family == "binomial"){

        rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "binomial",
                               cluster = cluster, cvControl = CV.control,
                               SL.library = list(c("SL.rpart", "All"),
                                                 c("SL.glm", "All"),
                                                 c("SL.gam", "All"),
                                                 c("SL.svm", "All"),
                                                 c("SL.glmnet", "All")))

        stopCluster(cluster)
        gc()
      }

    }else{

      if(family == "gaussian"){

        rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "gaussian",
                               cluster = cluster, cvControl = CV.control,
                               SL.library = library.SL)

        stopCluster(cluster)
        gc()

      }

      if(family == "binomial"){

        rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "binomial",
                               cluster = cluster, cvControl = CV.control,
                               SL.library = library.SL)

        stopCluster(cluster)
        gc()

      }

    }

  } else {

    if(is.null(library.SL)){

      if(family == "gaussian"){

        rc <- SuperLearner(X = x, Y = y, family = "gaussian",
                           cvControl = CV.control,
                           SL.library = list(c("SL.rpart", "All"),
                                             c("SL.glm", "All"),
                                             c("SL.gam", "All"),
                                             c("SL.svm", "All"),
                                             c("SL.glmnet", "All")))

      }

      if(family == "binomial"){

        rc <- SuperLearner(X = x, Y = y, family = "binomial",
                           cvControl = CV.control,
                           SL.library = list(c("SL.rpart", "All"),
                                             c("SL.glm", "All"),
                                             c("SL.gam", "All"),
                                             c("SL.svm", "All"),
                                             c("SL.glmnet", "All")))

      }

    }else{

      if(family == "gaussian"){

        rc <- SuperLearner(X = x, Y = y, family = "gaussian",
                           cvControl = CV.control,
                           SL.library = library.SL)

      }

      if(family == "binomial"){

        rc <- SuperLearner(X = x, Y = y, family = "binomial",
                           cvControl = CV.control,
                           SL.library = library.SL)

      }


    }

  }

  return(rc)

}

#' Predict final epidemic sizes using a trained SuperLearner
#'
#' @description
#' `SL_pred()` predicts the final epidemic sizes for locations given a trained SuperLearner objected (generated by `SL_fit()`) and covariates to be input into the SuperLearner.
#'
#'
#' @param mdl A `SuperLearner` object providing the trained ensemble model.
#' @param x A dataframe providing the variables/features for each location.
#'
#' @return A numeric vector of predicted epidemic final sizes for each location.
#'
#' @export
#'
#' @examples
#'
#' SL_pred(mdl = SL_model, x = my_vars)
#'
SL_pred <- function(mdl,
                    x) {

  pred <- pmax(0, predict(mdl, onlySL = T, newdata = x)$pred)

  return(pred)

}

#' Fit a generalized linear model
#'
#' @description
#' `linear_fit()` fits a generalized linear model based on location covariates and the observed cases at each location.
#'
#' @param x A data frame of covariates where each observation is a location
#' @param y A numeric vector object providing the outcome
#'
#' @return A glm object containing the model fit
#'
#' @export
#'
linear_fit <- function(x,
                       y,
                       stat.family = "gaussian") {

  tmp <- data.frame(y = y, x)

  rc <- glm(as.formula(paste0("y ~", paste(names(x), collapse = "+"))), data = tmp, family = stat.family)

  return(rc)

}

#' Predict epidemic size using a fitted glm model
#'
#' @description
#' `linear_pred()` predicts the final epidemic size for locations given a fitted glm model (generated using `linear_fit()`) and covariates to input into the model.
#'
#' @param mdl A glm object providing the model fit
#' @param x A data frame providing the covariates for each location
#'
#' @return A vector of predicted final sizes
#'
#' @export
#'
linear_pred <- function(mdl,
                        x) {

  pred <- pmax(0, predict(mdl, newdata = x))

  return(pred)

}
