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

test_that("plot_true_curves() isolated NA value in week", {

  epi_data_tmp <- epi_data
  epi_data_tmp$week[15] <- NA

  expect_s3_class(plot_true_curves(epi_data_tmp, week, district, cases), "ggplot")

})

test_that("plot_true_curves() isolated NA value in plot.group",{

  epi_data_tmp <- epi_data
  epi_data_tmp$district[4] <- NA

  expect_s3_class(plot_true_curves(epi_data_tmp, week, district, cases), "ggplot")

})

test_that("plot_true_curves() isolated NA value in count",{

  epi_data_tmp <- epi_data
  epi_data_tmp$cases[24] <- NA

  expect_s3_class(plot_true_curves(epi_data_tmp, week, district, cases), "ggplot")

})

test_that("plot_true_curves() empty data.frame", {

  empty_df <- NULL

  expect_error(plot_true_curves(empty_df, week, district, cases),
               "`dat` is empty")

  empty_df2 <- data.frame()

  expect_error(plot_true_curves(empty_df2, week, district, cases),
               "`dat` is empty")

})


test_that("plot_true_curves() plot.group variable all NAs", {

  epi_data_tmp <- epi_data
  epi_data_tmp$plot.group <- NA

  expect_s3_class(plot_true_curves(epi_data_tmp, week, district, cases),
                  "ggplot")

})

test_that("plot_true_curves() week variable all NAs", {

  epi_data_tmp <- epi_data
  epi_data_tmp$week <- NA

  expect_error(plot_true_curves(epi_data_tmp, week, district, cases),
               "`X` cannot contain only missing \\(NA\\) values")

})

test_that("plot_true_curves() count variable all NAs", {

  epi_data_tmp <- epi_data
  epi_data_tmp$cases <- NA

  expect_error(plot_true_curves(epi_data_tmp, week, district, cases),
               "`count` cannot contain only missing \\(NA\\) values")

})

test_that("plot_true_curves() missing dat argument", {

  expect_error(plot_true_curves(X = week, plot.group = district, count = cases),
               "`dat` is absent but must be supplied.")

})

test_that("plot_true_curves() missing X argument", {

  expect_error(plot_true_curves(dat = epi_data,
                                plot.group = district, count = cases),
               "`X` is absent but must be supplied.")

})

test_that("plot_true_curves() missing plot.group argument", {

  expect_error(plot_true_curves(dat = epi_data,
                                X = week, count = cases),
               "`plot.group` is absent but must be supplied.")

})

test_that("plot_true_curves() missing count argument", {

  expect_error(plot_true_curves(dat = epi_data, X = week,
                                plot.group = district),
               "`count` is absent but must be supplied.")

})

test_that("plot_true_curves() X variable not present in dat", {

  expect_error(plot_true_curves(dat = epi_data, X = year,
                                plot.group = district, count = cases),
               "`X` column\\(s\\) not found in data: \"year\"")

})

test_that("plot_true_curves() plot.group variable not present in dat", {

  expect_error(plot_true_curves(dat = epi_data, X = week,
                                plot.group = country, count = cases),
               "`plot.group` column\\(s\\) not found in data: \"country\"")

})

test_that("plot_true_curves() count variable not present in dat", {

  expect_error(plot_true_curves(dat = epi_data, X = week,
                                plot.group = district, count = n_inc),
               "`count` column\\(s\\) not found in data: \"n_inc\"")

})

test_that("plot_true_curves() character/factor count variable", {

  epi_data_tmp <- epi_data
  epi_data_tmp$cases <- as.character(epi_data_tmp$cases)

  expect_error(plot_true_curves(dat = epi_data_tmp, X = week,
                                plot.group = district, count = cases),
               "`count` must be of class <numeric>")

  epi_data_tmp2 <- epi_data
  epi_data_tmp2$cases <- as.factor(epi_data_tmp2$cases)

  expect_error(plot_true_curves(dat = epi_data_tmp2, X = week,
                                plot.group = district, count = cases),
               "`count` must be of class <numeric>")

})

test_that("plot_true_curves() dat is not a dataframe", {

  epi_data_ls <- list(district = rep(c("A1", "B1", "C1"), each = 10),
                      week = rep(1:10, times = 3),
                      cases = sample(1:15, size = 30, replace = T))

  expect_error(plot_true_curves(dat = epi_data_ls, X = week,
                                plot.group = district, count = cases),
               "`dat` must be of class <data.frame>")

})

test_that("plot_true_curves() multiple inputs for X variable", {

  epi_data_tmp <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                             week = rep(1:10, times = 3),
                             week2 = rep(1:10, times = 3),
                             cases = sample(1:15, size = 30, replace = T))

  expect_error(plot_true_curves(dat = epi_data_tmp, X = c(week, week2),
                                plot.group = district, count = cases),
               "`X` must be a single column name, not `c\\(week, week2\\)`")

})

test_that("plot_true_curves() multiple inputs for plot.groups variable", {

  epi_data_tmp <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                             country = rep(c("A", "B", "C"), each = 10),
                             week = rep(1:10, times = 3),
                             cases = sample(1:15, size = 30, replace = T))

  expect_error(plot_true_curves(dat = epi_data_tmp, X = week,
                                plot.group = c(district, country),
                                count = cases),
               "`plot.group` must be a single column name, not `c\\(district, country\\)`")

})

test_that("plot_true_curves() multiple inputs for count variable", {

  epi_data_tmp <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                             week = rep(1:10, times = 3),
                             cases = sample(1:15, size = 30, replace = T),
                             cases2 = sample(1:15, size = 30, replace = T))

  expect_error(plot_true_curves(dat = epi_data_tmp, X = week,
                                plot.group = district, count = c(cases, cases2)),
               "`count` must be a single column name, not `c\\(cases, cases2\\)`")

})

test_that("plot_true_curves() handles non-numeric classes for X", {

  epi_data_tmp <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                             week = as.Date(rep(1:10, times = 3)),
                             cases = sample(1:15, size = 30, replace = T))

  expect_no_error(plot_true_curves(dat = epi_data_tmp, X = week,
                                   plot.group = district, count = cases))

  epi_data_tmp2 <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                             week = as.factor(rep(1:10, times = 3)),
                             cases = sample(1:15, size = 30, replace = T))

  expect_no_error(plot_true_curves(dat = epi_data_tmp2, X = week,
                                   plot.group = district, count = cases))

  epi_data_tmp3 <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                             week = as.character(rep(1:10, times = 3)),
                             cases = sample(1:15, size = 30, replace = T))

  expect_no_error(plot_true_curves(dat = epi_data_tmp3, X = week,
                                   plot.group = district, count = cases))

})

test_that("plot_true_curves() handles non_numeric classes for plot.groups", {

  epi_data_tmp <- data.frame(district = as.factor(rep(c("A1", "B1", "C1"), each = 10)),
                             week = rep(1:10, times = 3),
                             cases = sample(1:15, size = 30, replace = T))

  expect_no_error(plot_true_curves(dat = epi_data_tmp, X = week,
                                   plot.group = district, count = cases))

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

test_that("plot_mask_curves() isolated NA value in week", {

  epi_data_tmp <- epi_data
  epi_data_tmp$week[15] <- NA

  expect_s3_class(plot_mask_curves(epi_data_tmp, mask_data, week, district, cases),
                  "ggplot")

  mask_data_tmp <- mask_data
  mask_data_tmp$week[4] <- NA

  expect_s3_class(plot_mask_curves(epi_data, mask_data_tmp, week, district, cases),
                  "ggplot")

})

test_that("plot_true_curves() isolated NA value in plot.group",{

  epi_data_tmp <- epi_data
  epi_data_tmp$district[4] <- NA

  expect_s3_class(plot_true_curves(epi_data_tmp, week, district, cases), "ggplot")

})

test_that("plot_true_curves() isolated NA value in count",{

  epi_data_tmp <- epi_data
  epi_data_tmp$cases[24] <- NA

  expect_s3_class(plot_true_curves(epi_data_tmp, week, district, cases), "ggplot")

})

test_that("plot_true_curves() empty data.frame", {

  empty_df <- NULL

  expect_error(plot_true_curves(empty_df, week, district, cases),
               "`dat` is empty")

  empty_df2 <- data.frame()

  expect_error(plot_true_curves(empty_df2, week, district, cases),
               "`dat` is empty")

})


test_that("plot_true_curves() plot.group variable all NAs", {

  epi_data_tmp <- epi_data
  epi_data_tmp$plot.group <- NA

  expect_s3_class(plot_true_curves(epi_data_tmp, week, district, cases),
                  "ggplot")

})

test_that("plot_true_curves() week variable all NAs", {

  epi_data_tmp <- epi_data
  epi_data_tmp$week <- NA

  expect_error(plot_true_curves(epi_data_tmp, week, district, cases),
               "`X` cannot contain only missing \\(NA\\) values")

})

test_that("plot_true_curves() count variable all NAs", {

  epi_data_tmp <- epi_data
  epi_data_tmp$cases <- NA

  expect_error(plot_true_curves(epi_data_tmp, week, district, cases),
               "`count` cannot contain only missing \\(NA\\) values")

})

test_that("plot_true_curves() missing dat argument", {

  expect_error(plot_true_curves(X = week, plot.group = district, count = cases),
               "`dat` is absent but must be supplied.")

})

test_that("plot_true_curves() missing X argument", {

  expect_error(plot_true_curves(dat = epi_data,
                                plot.group = district, count = cases),
               "`X` is absent but must be supplied.")

})

test_that("plot_true_curves() missing plot.group argument", {

  expect_error(plot_true_curves(dat = epi_data,
                                X = week, count = cases),
               "`plot.group` is absent but must be supplied.")

})

test_that("plot_true_curves() missing count argument", {

  expect_error(plot_true_curves(dat = epi_data, X = week,
                                plot.group = district),
               "`count` is absent but must be supplied.")

})

test_that("plot_true_curves() X variable not present in dat", {

  expect_error(plot_true_curves(dat = epi_data, X = year,
                                plot.group = district, count = cases),
               "`X` column\\(s\\) not found in data: \"year\"")

})

test_that("plot_true_curves() plot.group variable not present in dat", {

  expect_error(plot_true_curves(dat = epi_data, X = week,
                                plot.group = country, count = cases),
               "`plot.group` column\\(s\\) not found in data: \"country\"")

})

test_that("plot_true_curves() count variable not present in dat", {

  expect_error(plot_true_curves(dat = epi_data, X = week,
                                plot.group = district, count = n_inc),
               "`count` column\\(s\\) not found in data: \"n_inc\"")

})

test_that("plot_true_curves() character/factor count variable", {

  epi_data_tmp <- epi_data
  epi_data_tmp$cases <- as.character(epi_data_tmp$cases)

  expect_error(plot_true_curves(dat = epi_data_tmp, X = week,
                                plot.group = district, count = cases),
               "`count` must be of class <numeric>")

  epi_data_tmp2 <- epi_data
  epi_data_tmp2$cases <- as.factor(epi_data_tmp2$cases)

  expect_error(plot_true_curves(dat = epi_data_tmp2, X = week,
                                plot.group = district, count = cases),
               "`count` must be of class <numeric>")

})

test_that("plot_true_curves() dat is not a dataframe", {

  epi_data_ls <- list(district = rep(c("A1", "B1", "C1"), each = 10),
                      week = rep(1:10, times = 3),
                      cases = sample(1:15, size = 30, replace = T))

  expect_error(plot_true_curves(dat = epi_data_ls, X = week,
                                plot.group = district, count = cases),
               "`dat` must be of class <data.frame>")

})

test_that("plot_true_curves() multiple inputs for X variable", {

  epi_data_tmp <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                             week = rep(1:10, times = 3),
                             week2 = rep(1:10, times = 3),
                             cases = sample(1:15, size = 30, replace = T))

  expect_error(plot_true_curves(dat = epi_data_tmp, X = c(week, week2),
                                plot.group = district, count = cases),
               "`X` must be a single column name, not `c\\(week, week2\\)`")

})

test_that("plot_true_curves() multiple inputs for plot.groups variable", {

  epi_data_tmp <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                             country = rep(c("A", "B", "C"), each = 10),
                             week = rep(1:10, times = 3),
                             cases = sample(1:15, size = 30, replace = T))

  expect_error(plot_true_curves(dat = epi_data_tmp, X = week,
                                plot.group = c(district, country),
                                count = cases),
               "`plot.group` must be a single column name, not `c\\(district, country\\)`")

})

test_that("plot_true_curves() multiple inputs for count variable", {

  epi_data_tmp <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                             week = rep(1:10, times = 3),
                             cases = sample(1:15, size = 30, replace = T),
                             cases2 = sample(1:15, size = 30, replace = T))

  expect_error(plot_true_curves(dat = epi_data_tmp, X = week,
                                plot.group = district, count = c(cases, cases2)),
               "`count` must be a single column name, not `c\\(cases, cases2\\)`")

})

test_that("plot_true_curves() handles non-numeric classes for X", {

  epi_data_tmp <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                             week = as.Date(rep(1:10, times = 3)),
                             cases = sample(1:15, size = 30, replace = T))

  expect_no_error(plot_true_curves(dat = epi_data_tmp, X = week,
                                   plot.group = district, count = cases))

  epi_data_tmp2 <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                              week = as.factor(rep(1:10, times = 3)),
                              cases = sample(1:15, size = 30, replace = T))

  expect_no_error(plot_true_curves(dat = epi_data_tmp2, X = week,
                                   plot.group = district, count = cases))

  epi_data_tmp3 <- data.frame(district = rep(c("A1", "B1", "C1"), each = 10),
                              week = as.character(rep(1:10, times = 3)),
                              cases = sample(1:15, size = 30, replace = T))

  expect_no_error(plot_true_curves(dat = epi_data_tmp3, X = week,
                                   plot.group = district, count = cases))

})

test_that("plot_true_curves() handles non_numeric classes for plot.groups", {

  epi_data_tmp <- data.frame(district = as.factor(rep(c("A1", "B1", "C1"), each = 10)),
                             week = rep(1:10, times = 3),
                             cases = sample(1:15, size = 30, replace = T))

  expect_no_error(plot_true_curves(dat = epi_data_tmp, X = week,
                                   plot.group = district, count = cases))

})

#-------------------#
#### plot_bias() ####
#-------------------#


