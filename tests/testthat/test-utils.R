#### TEST UTILS.R ####

devtools::load_all()

#-------------------#
#### get_trueK() ####
#-------------------#

## Outputs

epi_data <- data.frame(place = c(rep("A1", times = 5),
                                 rep("B1", times = 10)),
                       cases_n = c(c(1:5), c(1:10)))

expt_vals <- c(sum(1:5), sum(1:10))

# correct class?
test_that("get_trueK() correct output class", {

  expect_s3_class(get_trueK(dat = epi_data,
                            case.col = cases_n,
                            groups = place),
                  "data.frame")
})

# correct dimensions?

test_that("get_trueK() correct output dimensions", {

  uniq_combs <- length(unique(epi_data$place))

  expect_equal(dim(get_trueK(dat = epi_data,
                             case.col = cases_n,
                             groups = place)),
               c(uniq_combs, 2))


})

# correct value?

test_that("get_trueK() correct output values", {

  expect_equal(get_trueK(dat = epi_data,
                         case.col = cases_n,
                         groups = place)$K,
               expt_vals)

})

## Inputs

# inputs w/ NA's

test_that("get_trueK() isolated NA case.col input", {

  na_data1 <- data.frame(place = c(rep("A1", times = 5),
                                   rep("B1", times = 10)),
                         cases_n = c(c(1:4, NA), c(1:3, NA, 5:10)))

  expect_false(any(is.na(get_trueK(na_data1,
                                   case.col = cases_n,
                                   groups = place)$K)))

})

test_that("get_trueK() isolated NA in groups input", {

  na_data2 <- data.frame(place = c(rep("A1", times = 4),
                                   NA,
                                   rep("B1", times = 10)),
                         cases_n = c(c(1:5), c(1:10)))

  expect_error(get_trueK(na_data2,
                         case.col = cases_n,
                         groups = place),
               "`place` cannot contain missing \\(NA\\) values")

})

# pass empty dataframe

test_that("get_trueK() empty df input", {

  empty_df <- NULL

  expect_error(get_trueK(empty_df,
                         case.col = cases_n,
                         groups = place),
               "`dat` is empty")

  empty_df2 <- data.frame()

  expect_error(get_trueK(empty_df2,
                        case.col = cases_n,
                        groups = place),
               "`dat` is empty")

})

# missing inputs

test_that("get_trueK() missing arguments", {

  expect_error(get_trueK(epi_data,
                         case.col = cases_n),
               "`groups` is absent but must be supplied.")

})

test_that("get_trueK() groups column(s) not present in dataframe", {

  expect_error(get_trueK(epi_data,
                         case.col = cases_n,
                         groups = "country"),
               "`groups` column\\(s\\) not found in data: \"country\"")

})

test_that("get_trueK() case.col column not present in dataframe", {

  expect_error(get_trueK(epi_data,
                         case.col = "my_cases",
                         groups = place),
               "`case.col` column\\(s\\) not found in data: \"my_cases\"")

})

# wrong class of input

test_that("get_trueK() character case.col input", {

  epi_data_tmp <- data.frame(place = c(rep("A1", times = 5),
                                       rep("B1", times = 10)),
                             cases_n = as.character(c(c(1:5), c(1:10))))

  expect_error(get_trueK(epi_data_tmp,
                         case.col = cases_n,
                         groups = place),
               "`case.col` must be of class <numeric>")

})

test_that("get_trueK() factor case.col input", {

  epi_data_tmp <- data.frame(place = c(rep("A1", times = 5),
                                       rep("B1", times = 10)),
                             cases_n = as.factor(c(c(1:5), c(1:10))))

  expect_error(get_trueK(epi_data_tmp,
                         case.col = cases_n,
                         groups = place),
               "`case.col` must be of class <numeric>")

})

# wrong type of input

test_that("get_trueK() dat is not a dataframe", {

  list_tmp <- list(place = c(rep("A1", times = 5),
                             rep("B1", times = 10)),
                   cases_n = c(c(1:5), c(1:10)))

  expect_error(get_trueK(list_tmp,
                         case.col = cases_n,
                         groups = place),
               "`dat` must be of class <data.frame>")

})
