#### TEST PLOT_FUNCS.R ####

devtools::load_all()

epi_data <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                       week = rep(1:10, times = 3),
                       cases = sample(1:15, size = 30, replace = T))

#--------------------------#
#### plot_true_curves() ####
#--------------------------#

test_that("plot_true_curves() returns a ggplot object", {

  expect_s3_class(plot_true_curves(epi_data, week, district, cases), "ggplot")
  expect_s3_class(plot_true_curves(epi_data, week, district, cases, legend = T),
                  "ggplot")

})

#--------------------------#
#### plot_mask_curves() ####
#--------------------------#

mask_data <- epi_data[c(1:5, 11:15, 21:25),]

test_that("plot_mask_curves() returns a ggplot object", {

  expect_s3_class(plot_mask_curves(epi_data, mask_data, week, district, cases),
                  "ggplot")
  expect_s3_class(plot_mask_curves(epi_data, mask_data, week, district, cases,
                                   legend = T),
                  "ggplot")

})

#-------------------#
#### plot_bias() ####
#-------------------#

K_vals = data.frame(loc = c("A1", "A2", "A3", "B1", "B2", "B3"),
                   K = c(15, 28, 78, 55, 21, 45))
iter_res = data.frame(A1 = c(25, 23, 20, 18, 17, 16),
                      A2 = c(40, 38, 37, 32, 31, 30),
                      A3 = c(50, 53, 61, 66, 72, 73),
                      B1 = c(35, 38, 39, 51, 52, 54),
                      B2 = c(40, 34, 33, 29, 27, 24),
                      B3 = c(62, 60, 59, 53, 49, 47))

test_that("plot_bias() returns a ggplot object", {

  expect_s3_class(plot_bias(trueK = K_vals, iter.results = iter_res),
                  "ggplot")

})
