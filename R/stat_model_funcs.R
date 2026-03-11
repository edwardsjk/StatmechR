## STAT MODEL FXNS ##

#' Train a SuperLearner ensemble in parallel
#'
#' @param x A dataframe providing the variables for each location.
#' @param y A numeric vector providing the outcome for each location.
#' @param cores A numeric object providing the number of cores to be used in the parallel process. The defaultt value is `1`.
#' @param family A character object providing the error distribution. Currenly only `gaussian` (the default) and `binomial` are supported.
#' @param library.SL A list object indicating the SL algorithms and screeners to be used in the Superlearner. Using the option `NULL` which will use the default set of algorithms and screeners.
#'
#' @return A `SuperLearner` object containing the trained model
#'
#' @importFrom SuperLearner snowSuperLearner
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterSetRNGStream
#' @importFrom parallel stopCluster
#'
#' @export
#'
#' @examples
#'
#' stat.mdl.sl.fit.para(x = my_vars, y = epi_size,
#'                      cores = 4, family = "gaussian")
#'
stat.mdl.sl.fit.para <- function(x, y, cores = 1, family = "gaussian", library.SL = NULL) {

  cluster = makeCluster(cores)
  clusterEvalQ(cluster,
               {library(SuperLearner)
                 options(mc.cores = 1)})
  clusterSetRNGStream(cluster, 1)

  if(is.null(library.SL)){

    if(family == "gaussian"){

      rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "gaussian",
                             cluster = cluster, cvControl = list(V = 10),
                             SL.library = list(c("SL.rpart", "screen.glmnet"),
                                               c("SL.glm", "screen.glmnet"),
                                               c('SL.gam', "screen.glmnet"),
                                               c("SL.glmnet", "screen.glmnet"),
                                               c("SL.svm", "screen.glmnet"),
                                               "SL.glmnet"))

      stopCluster(cluster)
      gc()

    }

    if(family == "binomial"){

      rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "binomial",
                             cluster = cluster, cvControl = list(V = 10),
                             SL.library = list(c("SL.rpart", "screen.glmnet"),
                                               c("SL.glm", "screen.glmnet"),
                                               c('SL.gam', "screen.glmnet"),
                                               c("SL.glmnet", "screen.glmnet"),
                                               c("SL.svm", "screen.glmnet"),
                                               "SL.glmnet"))

      stopCluster(cluster)
      gc()
    }

  }else{

    if(family == "gaussian"){

      rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "gaussian",
                             cluster = cluster, cvControl = list(V = 10),
                             SL.library = library.SL)

      stopCluster(cluster)
      gc()

    }

    if(family == "binomial"){

      rc <- snowSuperLearner(X = x, Y = y, newX = x, family = "binomial",
                             cluster = cluster, cvControl = list(V = 10),
                             SL.library = library.SL)

      stopCluster(cluster)
      gc()

    }

  }

  return(rc)

}

#' Train and fit a superlearner
#'
#' @param x A dataframe providing the variables for each location.
#' @param y A numeric vector providing the outcome for each location.
#' @param family A character object providing the error distribution. Currenly only `gaussian` (the default) and `poisson` are supported.
#' @param library.SL A list object indicating the SL algorithms and screeners to be used in the Superlearner. Using the option `NULL` which will use the default set of algorithms and screeners.
#'
#' @return A `SuperLearner` object containing the trained model
#'
#' @importFrom SuperLearner SuperLearner
#'
#' @export
#'
#' @examples
#'
#' stat.mdl.sl.fit(x = my_vars, y = epi_size, family = "gaussian")
#'
stat.mdl.sl.fit <- function(x, y, family = "gaussian", library.SL = NULL) {

  if(is.null(library.SL)){

    if(family == "gaussian"){

      rc <- SuperLearner(X = x, Y = y, family = "gaussian",
                         cvControl = list(V = 10),
                         SL.library = list(c("SL.rpart", "screen.glmnet"),
                                           c("SL.glm", "screen.glmnet"),
                                           c('SL.gam', "screen.glmnet"),
                                           c("SL.svm", "screen.glmnet"),
                                           "SL.glmnet"))

    }

    if(family == "binomial"){

      rc <- SuperLearner(X = x, Y = y, family = "binomial",
                         cvControl = list(V = 10),
                         SL.library = list(c("SL.rpart", "screen.glmnet"),
                                           c("SL.glm", "screen.glmnet"),
                                           c('SL.gam', "screen.glmnet"),
                                           c("SL.svm", "screen.glmnet"),
                                           "SL.glmnet"))

    }

  }else{

    if(family == "gaussian"){

      rc <- SuperLearner(X = x, Y = y, family = "gaussian",
                         cvControl = list(V = 10),
                         SL.library = library.SL)

    }

    if(family == "binomial"){

      rc <- SuperLearner(X = x, Y = y, family = "binomial",
                         cvControl = list(V = 10),
                         SL.library = library.SL)

    }


  }

  return(rc)

}

#' Predict final epidemic sizes using a trained SuperLearner
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
#' stat.mdl.sl.pred(mdl = SL_model, x = my_vars)
#'
stat.mdl.sl.pred <- function(mdl, x) {

  pred <- pmax(0, predict(mdl, onlySL = T, newdata = x)$pred)

  return(pred)

}

#' Function for fitting a simple linear statistical model appropriate for
#' passing in to em_func
#'
#' @param x the data frame of covariates to fit based on
#' @param y outcome
#' @return a fit model
#' @export
#'
stat.mdl.lin.fit <- function(x, y) {
  tmp <- data.frame(y = y, x)
  rc <- glm(as.formula(paste0("y ~", paste(names(x), collapse = "+"))), data = tmp, family = "gaussian")
  return(rc)
}

#'Function for doing the predict for doing the predict for the
#' results from stat.mdl.lin.fit
#'
#' @param mdl the fit model
#' @param x covariate matrix
#' @return a vector of predicted final sizes
#' @export
#'
stat.mdl.lin.pred <- function(mdl, x) {
  pred <- pmax(0, predict(mdl, newdata = x))
  return(pred)
}
