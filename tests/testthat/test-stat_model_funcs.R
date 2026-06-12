## TEST STAT_MODEL_FUNCS.R ##

devtools::load_all()

cars_outcome <- unname(unlist(mtcars |>
                                select(mpg)))

cars_vars <- mtcars |>
  select(-mpg)

rownames(cars_vars) <- 1:nrow(cars_vars)

##################
#### SL_fit() ####
##################

test_that("SL_fit returns correct output class", {

  expect_s3_class(SL_fit(x = cars_vars, y = cars_outcome, CV_folds = 3),
                  "SuperLearner")

})

test_that("SL_fit returns the same variables that are passed to it", {

  expect_equal(SL_fit(x = cars_vars, y = cars_outcome, CV_folds = 3)$varNames,
               colnames(cars_vars))

})

test_that("SL_fit returns the correct algorithms that are passed to it", {

  expect_equal(SL_fit(x = cars_vars, y = cars_outcome, CV_folds = 3,
                      library.SL = list(c("SL.rpart", "screen.corP"),
                                        c("SL.glm", "screen.corP"),
                                        c("SL.gam", "screen.corP"),
                                        c("SL.svm", "screen.corP")))$SL.library$library$predAlgorithm,
               c("SL.rpart", "SL.glm", "SL.gam", "SL.svm"))

})

###################
#### SL_pred() ####
###################

test_that("SL_pred() gives accurate results when provided data with a strong signal", {

  set.seed(42)
  n <- 200
  x_test <- data.frame(x1 = abs(rnorm(n)), x2 = abs(rnorm(n)))
  y_test <- 2 * x_test$x1 + 3 * x_test$x2  # perfect linear relationship, no noise

  mod <- SL_fit(x = x_test, y = y_test, CV_folds = 3)
  preds <- SL_pred(mdl = mod, x = x_test)

  expect_gt(cor(preds, y_test), 0.95)

})

test_that("SL_pred() beats a simple baseline", {

  mod <- SL_fit(x = cars_vars, y = cars_outcome, CV_folds = 3)
  preds <- SL_pred(mdl = mod, x = cars_vars)

  naive_pred <- mean(cars_outcome)  # simplest possible prediction
  model_rmse <- sqrt(mean((preds - cars_outcome)^2))
  naive_rmse <- sqrt(mean((naive_pred - cars_outcome)^2))

  expect_lt(model_rmse, naive_rmse)

})

test_that("SL_pred() returns a vector of the correct length", {

  mod <- SL_fit(x = cars_vars, y = cars_outcome, CV_folds = 3)
  preds <- SL_pred(mdl = mod, x = cars_vars)

  expect_equal(length(preds), nrow(cars_vars))

})

test_that("SL_pred() returns a vector where all values are non-negative", {

  mod <- SL_fit(x = cars_vars, y = cars_outcome, CV_folds = 3)
  preds <- SL_pred(mdl = mod, x = cars_vars)

  expect_true(all(preds >= 0))

})

######################
#### linear_fit() ####
######################

test_that("linear_fit() returns the correct class of output", {

  expect_s3_class(linear_fit(x = cars_vars, y = cars_outcome),
                  "glm")

})

test_that('linear_fit() returns correct number of variables in model object', {

  expect_equal(colnames(linear_fit(x = cars_vars, y = cars_outcome)$model)[-1],
         colnames(cars_vars))

})

#######################
#### linear_pred() ####
#######################

test_that("linear_pred() gives accurate results when provided data with a strong signal", {

  set.seed(42)
  n <- 200
  x_test <- data.frame(x1 = abs(rnorm(n)), x2 = abs(rnorm(n)))
  y_test <- 2 * x_test$x1 + 3 * x_test$x2  # perfect linear relationship, no noise

  mod <- linear_fit(x = x_test, y = y_test)
  preds <- linear_pred(mdl = mod, x = x_test)

  expect_gt(cor(preds, y_test), 0.95)

})

test_that("linear_pred() beats a simple baseline", {

  mod <- linear_fit(x = cars_vars, y = cars_outcome)
  preds <- linear_pred(mdl = mod, x = cars_vars)

  naive_pred <- mean(cars_outcome)  # simplest possible prediction
  model_rmse <- sqrt(mean((preds - cars_outcome)^2))
  naive_rmse <- sqrt(mean((naive_pred - cars_outcome)^2))

  expect_lt(model_rmse, naive_rmse)

})

test_that("linear()_pred returns a vector of the correct length", {

  mod <- linear_fit(x = cars_vars, y = cars_outcome)
  preds <- linear_pred(mdl = mod, x = cars_vars)

  expect_equal(length(preds), nrow(cars_vars))

})

test_that("linear_pred() returns a vector where all values are non-negative", {

  mod <- linear_fit(x = cars_vars, y = cars_outcome)
  preds <- linear_pred(mdl = mod, x = cars_vars)

  expect_true(all(preds >= 0))

})
