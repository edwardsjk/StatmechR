## TEST ERROR_FUNCS.R ##

devtools::load_all()

ec <- c(10, 10, 20, 30, 40, 50)
pred <- c(10, 20, 30, 50, 60, 50)
estK <- 300
prior <- 400

#-----------------------#
#### poisson_error() ####
#-----------------------#

test_that("poisson_error() gives correct error result", {

  err1 <- sum(dpois(ec, pred, log = TRUE))
  err2 <- dnorm(log10(estK), 0, 1, log = TRUE)

  expect_equal(poisson_error(ec, pred, estK),
               (err1 + err2))

})

test_that("poisson_error() handles 0 in ec", {

  ec2 <- c(0, 0, 0, 10, 20, 30)
  err1 <- sum(dpois(ec2, pred, log = TRUE))
  err2 <- dnorm(log10(estK), 0, 1, log = TRUE)

  expect_equal(poisson_error(ec2, pred, estK),
               (err1 + err2))

})

test_that("poisson_error() handles 0 in pred", {

  pred2 <- c(0, 0, 0, 10, 20, 30)
  err1 <- sum(dpois(ec, pred2, log = TRUE))
  err2 <- dnorm(log10(estK), 0, 1, log = TRUE)

  expect_equal(poisson_error(ec, pred2, estK),
               (err1 + err2))

})

test_that("poisson_error() handles 0 in estK", {

  err1 <- sum(dpois(ec, pred, log = TRUE))
  err2 <- dnorm(log10(0), 0, 1, log = TRUE)

  expect_equal(poisson_error(ec, pred, 0),
               (err1 + err2))

})

test_that("poisson_error() handles 0 in everything", {

  ec2 <- c(0, 0, 0, 0, 0, 0)
  pred2 <- c(0, 0, 0, 0, 0, 0)

  err1 <- sum(dpois(ec2, pred2, log = TRUE))
  err2 <- dnorm(log10(0), 0, 1, log = TRUE)

  expect_equal(poisson_error(ec2, pred2, 0),
               (err1 + err2))

})

#------------------------#
#### poisson_error2() ####
#------------------------#

test_that("poisson_error2() gives correct error result", {

  err1 <- sum(dpois(ec, pred, log = TRUE))
  err2 <- dnorm(log10(estK), log10(sum(ec)), 1, log = TRUE)

  expect_equal(poisson_error2(ec, pred, estK),
               (err1 + err2))

})

test_that("poisson_error2() handles some 0 in ec", {

  ec2 <- c(0, 0, 0, 10, 20, 30)
  err1 <- sum(dpois(ec2, pred, log = TRUE))
  err2 <- dnorm(log10(estK), log10(sum(ec2)), 1, log = TRUE)

  expect_equal(poisson_error2(ec2, pred, estK),
               (err1 + err2))

})

test_that("poisson_error2() handles some 0 in pred", {

  pred2 <- c(0, 0, 0, 10, 20, 30)
  err1 <- sum(dpois(ec, pred2, log = TRUE))
  err2 <- dnorm(log10(estK), log10(sum(ec)), 1, log = TRUE)

  expect_equal(poisson_error2(ec, pred2, estK),
               (err1 + err2))

})

test_that("poisson_error2() handles 0 in estK", {

  err1 <- sum(dpois(ec, pred, log = TRUE))
  err2 <- dnorm(log10(0), log10(sum(ec)), 1, log = TRUE)

  expect_equal(poisson_error2(ec, pred, 0),
               (err1 + err2))

})

test_that("poisson_error2() handles 0 in pred AND estK", {

  pred2 <- c(0, 0, 0, 0, 0, 0)

  err1 <- sum(dpois(ec, pred2, log = TRUE))
  err2 <- dnorm(log10(0), log10(sum(ec)), 1, log = TRUE)

  expect_equal(poisson_error2(ec, pred2, 0),
               (err1 + err2))

})

#------------------------#
#### gaussian_error() ####
#------------------------#

test_that("gaussian_error() gives correct error result", {

  err1 <- sum(dnorm(ec, pred, 1, log = TRUE))
  err2 <- dnorm(log10(estK), 0, 1, log = TRUE)

  expect_equal(gaussian_error(ec, pred, estK),
               (err1 + err2))

})

test_that("gaussian_error() handles 0 in ec", {

  ec2 <- c(0, 0, 0, 10, 20, 30)
  err1 <- sum(dnorm(ec2, pred, 1, log = TRUE))
  err2 <- dnorm(log10(estK), 0, 1, log = TRUE)

  expect_equal(gaussian_error(ec2, pred, estK),
               (err1 + err2))

})

test_that("gaussian_error() handles 0 in pred", {

  pred2 <- c(0, 0, 0, 10, 20, 30)
  err1 <- sum(dnorm(ec, pred2, 1, log = TRUE))
  err2 <- dnorm(log10(estK), 0, 1, log = TRUE)

  expect_equal(gaussian_error(ec, pred2, estK),
               (err1 + err2))

})

test_that("gaussian_error() handles 0 in estK", {

  err1 <- sum(dnorm(ec, pred, 1, log = TRUE))
  err2 <- dnorm(log10(0), 0, 1, log = TRUE)

  expect_equal(gaussian_error(ec, pred, 0),
               (err1 + err2))

})

test_that("gaussian_error() handles 0 in everything", {

  pred2 <- c(0, 0, 0, 0, 0, 0)

  err1 <- sum(dnorm(ec, pred2, 1, log = TRUE))
  err2 <- dnorm(log10(0), 0, 1, log = TRUE)

  expect_equal(gaussian_error(ec, pred2, 0),
               (err1 + err2))

})

#-------------------------#
#### gaussian_error2() ####
#-------------------------#

test_that("gaussian_error2() gives correct error result", {

  err1 <- sum(dnorm(ec, pred, 1, log = TRUE))
  err2 <- dnorm(log10(estK), log10(sum(ec)), 1, log = TRUE)

  expect_equal(gaussian_error2(ec, pred, estK),
               (err1 + err2))

})

test_that("gaussian_error2() handles 0 in ec", {

  ec2 <- c(0, 0, 0, 10, 20, 30)
  err1 <- sum(dnorm(ec2, pred, 1, log = TRUE))
  err2 <- dnorm(log10(estK), log10(sum(ec2)), 1, log = TRUE)

  expect_equal(gaussian_error2(ec2, pred, estK),
               (err1 + err2))

})

test_that("gaussian_error2() handles 0 in pred", {

  pred2 <- c(0, 0, 0, 10, 20, 30)
  err1 <- sum(dnorm(ec, pred2, 1, log = TRUE))
  err2 <- dnorm(log10(estK), log10(sum(ec)), 1, log = TRUE)

  expect_equal(gaussian_error2(ec, pred2, estK),
               (err1 + err2))

})

test_that("gaussian_error2() handles 0 in estK", {

  err1 <- sum(dnorm(ec, pred, 1, log = TRUE))
  err2 <- dnorm(log10(0), log10(sum(ec)), 1, log = TRUE)

  expect_equal(gaussian_error(ec, pred, 0),
               (err1 + err2))

})

test_that("gaussian_error2() handles 0 in everything", {

  pred2 <- c(0, 0, 0, 0, 0, 0)

  err1 <- sum(dnorm(ec, pred2, 1, log = TRUE))
  err2 <- dnorm(log10(0), log10(sum(ec)), 1, log = TRUE)

  expect_equal(gaussian_error(ec, pred2, 0),
               (err1 + err2))

})

#-----------------#
#### poispen() ####
#-----------------#

test_that("poispen() gives correct error result", {

  err <- dpois(round(estK), round(prior), log = TRUE)

  expect_equal(poispen(estK, prior),
               err)

})

test_that("poispen() handles 0 estK", {

  err <- dpois(round(0), round(prior), log = TRUE)

  expect_equal(poispen(0, prior),
               err)

})

test_that("poispen() handles 0 prior", {

  err <- dpois(round(estK), round(0), log = TRUE)

  expect_equal(poispen(estK, 0),
               err)

})

test_that("poispen() handles 0 estK AND 0 prior", {

  err <- dpois(round(0), round(0), log = TRUE)

  expect_equal(poispen(0, 0),
               err)

})

#-----------------#
#### sqrtpen() ####
#-----------------#

test_that("sqrtpen() gives correct error value", {

  err <- dnorm(sqrt(abs((estK) - prior)), 0, 1, log = TRUE)

  expect_equal(sqrtpen(estK, prior),
               err)

})

test_that("sqrtpen() handles 0 estK", {

  err <- dnorm(sqrt(abs((0) - prior)), 0, 1, log = TRUE)

  expect_equal(sqrtpen(0, prior),
               err)

})

test_that("sqrtpen() handles 0 prior", {

  err <- dnorm(sqrt(abs((estK) - 0)), 0, 1, log = TRUE)

  expect_equal(sqrtpen(estK, 0),
               err)

})

test_that("sqrtpen() handles 0 estK AND 0 prior", {

  err <- dnorm(sqrt(abs((0) - 0)), 0, 1, log = TRUE)

  expect_equal(sqrtpen(0, 0),
               err)

})

#----------------------#
#### sqrtpen_diff() ####
#----------------------#

test_that("sqrtpen_diff() gives correct error value", {

  err <- dnorm(sqrt(abs((estK) - prior)), 0, 2, log = TRUE)

  expect_equal(sqrtpen_diff(estK, prior),
               err)

})

test_that("sqrtpen_diff() handles 0 estK", {

  err <- dnorm(sqrt(abs((0) - prior)), 0, 2, log = TRUE)

  expect_equal(sqrtpen_diff(0, prior),
               err)

})

test_that("sqrtpen_diff() handles 0 prior", {

  err <- dnorm(sqrt(abs((estK) - 0)), 0, 2, log = TRUE)

  expect_equal(sqrtpen_diff(estK, 0),
               err)

})

test_that("sqrtpen_diff() handles 0 estK AND 0 prior", {

  err <- dnorm(sqrt(abs((0) - 0)), 0, 2, log = TRUE)

  expect_equal(sqrtpen_diff(0, 0),
               err)

})

#------------------#
#### gausspen() ####
#------------------#

test_that("gausspen() gives correct error value", {

  err <- dnorm((abs((estK) - prior)), 0, 1, log = TRUE)

  expect_equal(gausspen(estK, prior),
               err)

})

test_that("gausspen() handles 0 estK", {

  err <- dnorm((abs((0) - prior)), 0, 1, log = TRUE)

  expect_equal(gausspen(0, prior),
               err)

})

test_that("gausspen() handles 0 prior", {

  err <- dnorm((abs((estK) - 0)), 0, 1, log = TRUE)

  expect_equal(gausspen(estK, 0),
               err)

})

test_that("gausspen() handles 0 estK and 0 prior", {

  err <- dnorm((abs((0) - 0)), 0, 1, log = TRUE)

  expect_equal(gausspen(0, 0),
               err)

})
