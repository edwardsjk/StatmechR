## TEST ERROR_FUNCS.R ##

devtools::load_all()

ec <- c(10, 10, 20, 30, 40, 50)
pred <- c(10, 20, 30, 50, 60, 50)

#-----------------------#
#### poisson_error() ####
#-----------------------#

test_that("poisson_error() correct output class", {

  expect_equal(class(poisson_error(ec, pred, sum(pred))),
                  "numeric")

})

test_that("poisson_error() correct shape/length of output", {

  expect_vector(poisson_error(ec, pred, sum(pred)), size = 1)

})

test_that("poisson_error() ec contains NA", {

  ec_tmp <- c(10, 10, NA, 30, 40, 50)

  expect_error(poisson_error(ec_tmp, pred, sum(pred)), "`ec` cannot contain missing \\(NA\\) values")

})

test_that("poisson_error() pred contains NA", {

  pred_tmp <- c(10, 20, 30, 50, NA, 50)

  expect_error(poisson_error(ec, pred_tmp, sum(pred)), "`pred` cannot contain missing \\(NA\\) values")

})

test_that("poisson_error() estK contains NA", {

  estK_tmp <- NA

  expect_error(poisson_error(ec, pred, estK_tmp), "`estK` cannot contain missing \\(NA\\) values")

})

test_that("poisson_error() empty ec input", {

  ec_tmp <- NULL

  expect_error(poisson_error(ec_tmp, pred, sum(pred)), "`ec` is empty")

})

test_that("poisson_error() empty pred input", {

  pred_tmp <- NULL

  expect_error(poisson_error(ec, pred_tmp, sum(pred)), "`pred` is empty")

})

test_that("poisson_error() empty estK input", {

  estK_tmp <- NULL

  expect_error(poisson_error(ec, pred, estK_tmp), "`estK` is empty")

})

test_that("poisson_error() character/factor ec input", {

  ec1 <- as.character(ec)
  ec2 <- as.factor(ec)

  expect_error(poisson_error(ec1, pred, sum(pred)), "`ec` must be of class <numeric>")
  expect_error(poisson_error(ec2, pred, sum(pred)), "`ec` must be of class <numeric>")

})

test_that("poisson_error() character/factor pred input", {

  pred1 <- as.character(pred)
  pred2 <- as.factor(pred)

  expect_error(poisson_error(ec, pred1, sum(pred)), "`pred` must be of class <numeric>")
  expect_error(poisson_error(ec, pred2, sum(pred)), "`pred` must be of class <numeric>")

})

test_that("poisson_error() character/factor estk input", {

  estK1 <- as.character(sum(pred))
  estK2 <- as.factor(sum(pred))

  expect_error(poisson_error(ec, pred, estK1), "`estK` must be of class <numeric>")
  expect_error(poisson_error(ec, pred, estK2), "`estK` must be of class <numeric>")

})

test_that("poisson_error() non-positive estK input", {

  estK_tmp <- -100

  expect_error(poisson_error(ec, pred, estK_tmp),
               "`estK` must equal 0 or a positive value")

})

test_that("poisson_error() missing ec argument", {

  expect_error(poisson_error(pred = pred, estK = sum(pred)),
               "`ec` is absent but must be supplied.")

})

test_that("poisson_error() missing pred argument", {

  expect_error(poisson_error(ec = ec, estK = sum(pred)),
               "`pred` is absent but must be supplied.")

})

test_that("poisson_error() missing estK argument", {

  expect_error(poisson_error(ec = ec, pred = pred),
               "`estK` is absent but must be supplied.")

})

#------------------------#
#### poisson_error2() ####
#------------------------#

test_that("poisson_error2() correct output class", {

  expect_equal(class(poisson_error2(ec, pred, sum(pred))),
               "numeric")

})

test_that("poisson_error2() correct shape/length of output", {

  expect_vector(poisson_error2(ec, pred, sum(pred)), size = 1)

})

test_that("poisson_error2() ec contains NA", {

  ec_tmp <- c(10, 10, NA, 30, 40, 50)

  expect_error(poisson_error2(ec_tmp, pred, sum(pred)), "`ec` cannot contain missing \\(NA\\) values")

})

test_that("poisson_error2() pred contains NA", {

  pred_tmp <- c(10, 20, 30, 50, NA, 50)

  expect_error(poisson_error2(ec, pred_tmp, sum(pred)), "`pred` cannot contain missing \\(NA\\) values")

})

test_that("poisson_error2() estK contains NA", {

  estK_tmp <- NA

  expect_error(poisson_error2(ec, pred, estK_tmp), "`estK` cannot contain missing \\(NA\\) values")

})

test_that("poisson_error2() empty ec input", {

  ec_tmp <- NULL

  expect_error(poisson_error2(ec_tmp, pred, sum(pred)), "`ec` is empty")

})

test_that("poisson_error2() empty pred input", {

  pred_tmp <- NULL

  expect_error(poisson_error2(ec, pred_tmp, sum(pred)), "`pred` is empty")

})

test_that("poisson_error2() empty estK input", {

  estK_tmp <- NULL

  expect_error(poisson_error2(ec, pred, estK_tmp), "`estK` is empty")

})

test_that("poisson_error2() character/factor ec input", {

  ec1 <- as.character(ec)
  ec2 <- as.factor(ec)

  expect_error(poisson_error2(ec1, pred, sum(pred)), "`ec` must be of class <numeric>")
  expect_error(poisson_error2(ec2, pred, sum(pred)), "`ec` must be of class <numeric>")

})

test_that("poisson_error2() character/factor pred input", {

  pred1 <- as.character(pred)
  pred2 <- as.factor(pred)

  expect_error(poisson_error2(ec, pred1, sum(pred)), "`pred` must be of class <numeric>")
  expect_error(poisson_error2(ec, pred2, sum(pred)), "`pred` must be of class <numeric>")

})

test_that("poisson_error2() character/factor estk input", {

  estK1 <- as.character(sum(pred))
  estK2 <- as.factor(sum(pred))

  expect_error(poisson_error2(ec, pred, estK1), "`estK` must be of class <numeric>")
  expect_error(poisson_error2(ec, pred, estK2), "`estK` must be of class <numeric>")

})

test_that("poisson_error2() non-positive estK input", {

  estK_tmp <- -100

  expect_error(poisson_error2(ec, pred, estK_tmp), "`estK` must equal 0 or a positive value")

})

test_that("poisson_error2() missing ec argument", {

  expect_error(poisson_error2(pred = pred, estK = sum(pred)),
               "`ec` is absent but must be supplied.")

})

test_that("poisson_error2() missing pred argument", {

  expect_error(poisson_error2(ec = ec, estK = sum(pred)),
               "`pred` is absent but must be supplied.")

})

test_that("poisson_error2() missing estK argument", {

  expect_error(poisson_error2(ec = ec, pred = pred),
               "`estK` is absent but must be supplied.")

})

#------------------------#
#### gaussian_error() ####
#------------------------#

test_that("gaussian_error() correct output class", {

  expect_equal(class(gaussian_error(ec, pred, sum(pred))),
               "numeric")

})

test_that("gaussian_error() correct shape/length of output", {

  expect_vector(gaussian_error(ec, pred, sum(pred)), size = 1)

})

test_that("gaussian_error() ec contains NA", {

  ec_tmp <- c(10, 10, NA, 30, 40, 50)

  expect_error(gaussian_error(ec_tmp, pred, sum(pred)), "`ec` cannot contain missing \\(NA\\) values")

})

test_that("gaussian_error() pred contains NA", {

  pred_tmp <- c(10, 20, 30, 50, NA, 50)

  expect_error(gaussian_error(ec, pred_tmp, sum(pred)), "`pred` cannot contain missing \\(NA\\) values")

})

test_that("gaussian_error() estK contains NA", {

  estK_tmp <- NA

  expect_error(gaussian_error(ec, pred, estK_tmp), "`estK` cannot contain missing \\(NA\\) values")

})

test_that("gaussian_error() empty ec input", {

  ec_tmp <- NULL

  expect_error(gaussian_error(ec_tmp, pred, sum(pred)), "`ec` is empty")

})

test_that("gaussian_error() empty pred input", {

  pred_tmp <- NULL

  expect_error(gaussian_error(ec, pred_tmp, sum(pred)), "`pred` is empty")

})

test_that("gaussian_error() empty estK input", {

  estK_tmp <- NULL

  expect_error(gaussian_error(ec, pred, estK_tmp), "`estK` is empty")

})

test_that("gaussian_error() character/factor ec input", {

  ec1 <- as.character(ec)
  ec2 <- as.factor(ec)

  expect_error(gaussian_error(ec1, pred, sum(pred)), "`ec` must be of class <numeric>")
  expect_error(gaussian_error(ec2, pred, sum(pred)), "`ec` must be of class <numeric>")

})

test_that("gaussian_error() character/factor pred input", {

  pred1 <- as.character(pred)
  pred2 <- as.factor(pred)

  expect_error(gaussian_error(ec, pred1, sum(pred)), "`pred` must be of class <numeric>")
  expect_error(gaussian_error(ec, pred2, sum(pred)), "`pred` must be of class <numeric>")

})

test_that("gaussian_error() character/factor estk input", {

  estK1 <- as.character(sum(pred))
  estK2 <- as.factor(sum(pred))

  expect_error(gaussian_error(ec, pred, estK1), "`estK` must be of class <numeric>")
  expect_error(gaussian_error(ec, pred, estK2), "`estK` must be of class <numeric>")

})

test_that("gaussian_error() non-positive estK input", {

  estK_tmp <- -100

  expect_error(gaussian_error(ec, pred, estK_tmp), "`estK` must equal 0 or a positive value")

})

test_that("gaussian_error() missing ec argument", {

  expect_error(gaussian_error(pred = pred, estK = sum(pred)),
               "`ec` is absent but must be supplied.")

})

test_that("gaussian_error() missing pred argument", {

  expect_error(gaussian_error(ec = ec, estK = sum(pred)),
               "`pred` is absent but must be supplied.")

})

test_that("gaussian_error() missing estK argument", {

  expect_error(gaussian_error(ec = ec, pred = pred),
               "`estK` is absent but must be supplied.")

})

#-------------------------#
#### gaussian_error2() ####
#-------------------------#

test_that("gaussian_error2() correct output class", {

  expect_equal(class(gaussian_error2(ec, pred, sum(pred))),
               "numeric")

})

test_that("gaussian_error2() correct shape/length of output", {

  expect_vector(gaussian_error2(ec, pred, sum(pred)), size = 1)

})

test_that("gaussian_error2() ec contains NA", {

  ec_tmp <- c(10, 10, NA, 30, 40, 50)

  expect_error(gaussian_error2(ec_tmp, pred, sum(pred)), "`ec` cannot contain missing \\(NA\\) values")

})

test_that("gaussian_error2() pred contains NA", {

  pred_tmp <- c(10, 20, 30, 50, NA, 50)

  expect_error(gaussian_error2(ec, pred_tmp, sum(pred)), "`pred` cannot contain missing \\(NA\\) values")

})

test_that("gaussian_error2() estK contains NA", {

  estK_tmp <- NA

  expect_error(gaussian_error2(ec, pred, estK_tmp), "`estK` cannot contain missing \\(NA\\) values")

})

test_that("gaussian_error2() empty ec input", {

  ec_tmp <- NULL

  expect_error(gaussian_error2(ec_tmp, pred, sum(pred)), "`ec` is empty")

})

test_that("gaussian_error2() empty pred input", {

  pred_tmp <- NULL

  expect_error(gaussian_error2(ec, pred_tmp, sum(pred)), "`pred` is empty")

})

test_that("gaussian_error2() empty estK input", {

  estK_tmp <- NULL

  expect_error(gaussian_error2(ec, pred, estK_tmp), "`estK` is empty")

})

test_that("gaussian_error2() character/factor ec input", {

  ec1 <- as.character(ec)
  ec2 <- as.factor(ec)

  expect_error(gaussian_error2(ec1, pred, sum(pred)), "`ec` must be of class <numeric>")
  expect_error(gaussian_error2(ec2, pred, sum(pred)), "`ec` must be of class <numeric>")

})

test_that("gaussian_error2() character/factor pred input", {

  pred1 <- as.character(pred)
  pred2 <- as.factor(pred)

  expect_error(gaussian_error2(ec, pred1, sum(pred)), "`pred` must be of class <numeric>")
  expect_error(gaussian_error2(ec, pred2, sum(pred)), "`pred` must be of class <numeric>")

})

test_that("gaussian_error2() character/factor estk input", {

  estK1 <- as.character(sum(pred))
  estK2 <- as.factor(sum(pred))

  expect_error(gaussian_error2(ec, pred, estK1), "`estK` must be of class <numeric>")
  expect_error(gaussian_error2(ec, pred, estK2), "`estK` must be of class <numeric>")

})

test_that("gaussian_error2() non-positive estK input", {

  estK_tmp <- -100

  expect_error(gaussian_error2(ec, pred, estK_tmp), "`estK` must equal 0 or a positive value")

})

test_that("gaussian_error2() missing ec argument", {

  expect_error(gaussian_error2(pred = pred, estK = sum(pred)),
               "`ec` is absent but must be supplied.")

})

test_that("gaussian_error2() missing pred argument", {

  expect_error(gaussian_error2(ec = ec, estK = sum(pred)),
               "`pred` is absent but must be supplied.")

})

test_that("gaussian_error2() missing estK argument", {

  expect_error(gaussian_error2(ec = ec, pred = pred),
               "`estK` is absent but must be supplied.")

})

#-----------------#
#### poispen() ####
#-----------------#

prior <- 1000
estK <- 900

test_that("poispen() correct output class", {

  expect_equal(class(poispen(estK, prior)),
               "numeric")

})

test_that("poispen() correct shape/length of output", {

  expect_vector(poispen(estK, prior), size = 1)

})

test_that("poispen() estK contains NA", {

   estK_tmp <- NA

  expect_error(poispen(estK_tmp, prior), "`estK` cannot contain missing \\(NA\\) values")

})

test_that("poispen() prior contains NA", {

  prior_tmp <- NA

  expect_error(poispen(estK, prior_tmp), "`prior` cannot contain missing \\(NA\\) values")

})

test_that("poispen() empty estK input", {

  estK_tmp <- NULL

  expect_error(poispen(estK_tmp, prior), "`estK` is empty")

})

test_that("poispen() empty prior input", {

  prior_tmp <- NULL

  expect_error(poispen(estK, prior_tmp), "`prior` is empty")

})

test_that("poispen() character/factor estK input", {

  estK1 <- as.character(estK)
  estK2 <- as.factor(estK)

  expect_error(poispen(estK1, prior), "`estK` must be of class <numeric>")
  expect_error(poispen(estK2, prior), "`estK` must be of class <numeric>")

})

test_that("poispen() character/factor prior input", {

  prior1 <- as.character(prior)
  prior2 <- as.factor(prior)

  expect_error(poispen(estK, prior1), "`prior` must be of class <numeric>")
  expect_error(poispen(estK, prior2), "`prior` must be of class <numeric>")

})

test_that("poispen() non-positive estK input", {

  estK_tmp <- -100

  expect_error(poispen(estK_tmp, prior), "`estK` must equal 0 or a positive value")

})

test_that("poispen() non-positive prior input", {

  prior_tmp <- -100

  expect_error(poispen(estK, prior_tmp), "`prior` must equal 0 or a positive value")

})

test_that("poispen() missing estK argument", {

  expect_error(poispen(prior = prior),
               "`estK` is absent but must be supplied.")

})

test_that("poispen() missing prior argument", {

  expect_error(poispen(estK = estK),
               "`prior` is absent but must be supplied.")

})

#-----------------#
#### sqrtpen() ####
#-----------------#

test_that("sqrtpen() correct output class", {

  expect_equal(class(sqrtpen(estK, prior)),
               "numeric")

})

test_that("sqrtpen() correct shape/length of output", {

  expect_vector(sqrtpen(estK, prior), size = 1)

})

test_that("sqrtpen() estK contains NA", {

  estK_tmp <- NA

  expect_error(sqrtpen(estK_tmp, prior), "`estK` cannot contain missing \\(NA\\) values")

})

test_that("sqrtpen() prior contains NA", {

  prior_tmp <- NA

  expect_error(sqrtpen(estK, prior_tmp), "`prior` cannot contain missing \\(NA\\) values")

})

test_that("sqrtpen() empty estK input", {

  estK_tmp <- NULL

  expect_error(sqrtpen(estK_tmp, prior), "`estK` is empty")

})

test_that("sqrtpen() empty prior input", {

  prior_tmp <- NULL

  expect_error(sqrtpen(estK, prior_tmp), "`prior` is empty")

})

test_that("sqrtpen() character/factor estK input", {

  estK1 <- as.character(estK)
  estK2 <- as.factor(estK)

  expect_error(sqrtpen(estK1, prior), "`estK` must be of class <numeric>")
  expect_error(sqrtpen(estK2, prior), "`estK` must be of class <numeric>")

})

test_that("sqrtpen() character/factor prior input", {

  prior1 <- as.character(prior)
  prior2 <- as.factor(prior)

  expect_error(sqrtpen(estK, prior1), "`prior` must be of class <numeric>")
  expect_error(sqrtpen(estK, prior2), "`prior` must be of class <numeric>")

})

test_that("sqrtpen() non-positive estK input", {

  estK_tmp <- -100

  expect_error(sqrtpen(estK_tmp, prior), "`estK` must equal 0 or a positive value")

})

test_that("sqrtpen() non-positive prior input", {

  prior_tmp <- -100

  expect_error(sqrtpen(estK, prior_tmp), "`prior` must equal 0 or a positive value")

})

test_that("sqrtpen() missing estK argument", {

  expect_error(sqrtpen(prior = prior),
               "`estK` is absent but must be supplied.")

})

test_that("sqrtpen() missing prior argument", {

  expect_error(sqrtpen(estK = estK),
               "`prior` is absent but must be supplied.")

})

#----------------------#
#### sqrtpen_diff() ####
#----------------------#

test_that("sqrtpen_diff() correct output class", {

  expect_equal(class(sqrtpen_diff(estK, prior)),
               "numeric")

})

test_that("sqrtpen_diff() correct shape/length of output", {

  expect_vector(sqrtpen_diff(estK, prior), size = 1)

})

test_that("sqrtpen_diff() estK contains NA", {

  estK_tmp <- NA

  expect_error(sqrtpen_diff(estK_tmp, prior), "`estK` cannot contain missing \\(NA\\) values")

})

test_that("sqrtpen_diff() prior contains NA", {

  prior_tmp <- NA

  expect_error(sqrtpen_diff(estK, prior_tmp), "`prior` cannot contain missing \\(NA\\) values")

})

test_that("sqrtpen_diff() empty estK input", {

  estK_tmp <- NULL

  expect_error(sqrtpen_diff(estK_tmp, prior), "`estK` is empty")

})

test_that("sqrtpen_diff() empty prior input", {

  prior_tmp <- NULL

  expect_error(sqrtpen_diff(estK, prior_tmp), "`prior` is empty")

})

test_that("sqrtpen_diff() character/factor estK input", {

  estK1 <- as.character(estK)
  estK2 <- as.factor(estK)

  expect_error(sqrtpen_diff(estK1, prior), "`estK` must be of class <numeric>")
  expect_error(sqrtpen_diff(estK2, prior), "`estK` must be of class <numeric>")

})

test_that("sqrtpen_diff() character/factor prior input", {

  prior1 <- as.character(prior)
  prior2 <- as.factor(prior)

  expect_error(sqrtpen_diff(estK, prior1), "`prior` must be of class <numeric>")
  expect_error(sqrtpen_diff(estK, prior2), "`prior` must be of class <numeric>")

})

test_that("sqrtpen_diff() non-positive estK input", {

  estK_tmp <- -100

  expect_error(sqrtpen_diff(estK_tmp, prior), "`estK` must equal 0 or a positive value")

})

test_that("sqrtpen_diff() non-positive prior input", {

  prior_tmp <- -100

  expect_error(sqrtpen_diff(estK, prior_tmp), "`prior` must equal 0 or a positive value")

})

test_that("sqrtpen_diff() missing estK argument", {

  expect_error(sqrtpen_diff(prior = prior),
               "`estK` is absent but must be supplied.")

})

test_that("sqrtpen_diff() missing prior argument", {

  expect_error(sqrtpen_diff(estK = estK),
               "`prior` is absent but must be supplied.")

})

#------------------#
#### gausspen() ####
#------------------#

test_that("gausspen() correct output class", {

  expect_equal(class(gausspen(estK, prior)),
               "numeric")

})

test_that("gausspen() correct shape/length of output", {

  expect_vector(gausspen(estK, prior), size = 1)

})

test_that("gausspen() estK contains NA", {

  estK_tmp <- NA

  expect_error(gausspen(estK_tmp, prior), "`estK` cannot contain missing \\(NA\\) values")

})

test_that("gausspen() prior contains NA", {

  prior_tmp <- NA

  expect_error(gausspen(estK, prior_tmp), "`prior` cannot contain missing \\(NA\\) values")

})

test_that("gausspen() empty estK input", {

  estK_tmp <- NULL

  expect_error(gausspen(estK_tmp, prior), "`estK` is empty")

})

test_that("gausspen() empty prior input", {

  prior_tmp <- NULL

  expect_error(gausspen(estK, prior_tmp), "`prior` is empty")

})

test_that("gausspen() character/factor estK input", {

  estK1 <- as.character(estK)
  estK2 <- as.factor(estK)

  expect_error(gausspen(estK1, prior), "`estK` must be of class <numeric>")
  expect_error(gausspen(estK2, prior), "`estK` must be of class <numeric>")

})

test_that("gausspen() character/factor prior input", {

  prior1 <- as.character(prior)
  prior2 <- as.factor(prior)

  expect_error(gausspen(estK, prior1), "`prior` must be of class <numeric>")
  expect_error(gausspen(estK, prior2), "`prior` must be of class <numeric>")

})

test_that("gausspen() non-positive estK input", {

  estK_tmp <- -100

  expect_error(gausspen(estK_tmp, prior), "`estK` must equal 0 or a positive value")

})

test_that("gausspen() non-positive prior input", {

  prior_tmp <- -100

  expect_error(gausspen(estK, prior_tmp), "`prior` must equal 0 or a positive value")

})

test_that("gausspen() missing estK argument", {

  expect_error(gausspen(prior = prior),
               "`estK` is absent but must be supplied.")

})

test_that("gausspen() missing prior argument", {

  expect_error(gausspen(estK = estK),
               "`prior` is absent but must be supplied.")

})

