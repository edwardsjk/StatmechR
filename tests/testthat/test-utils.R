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


test_that("get_trueK() correct output values (1 group variables)", {

  expect_equal(get_trueK(dat = epi_data,
                         case.col = cases_n,
                         groups = place)$K,
               expt_vals)

})

epi_data2 <- data.frame(place = c(rep("A1", times = 5),
                                 rep("A2", times = 7),
                                 rep("A3", times = 12),
                                 rep("B1", times = 10),
                                 rep("B2", times = 6),
                                 rep("B3", times = 9)),
                       country = c(rep("Country A", times = 24),
                                   rep("Country B", times = 25)),
                       cases_n = c(1:5, 1:7, 1:12, 1:10, 1:6, 1:9))

expt_vals2 <- c(sum(1:5), sum(1:7), sum(1:12), sum(1:10), sum(1:6), sum(1:9))

test_that("get_trueK() correct output values (2 group variables)", {

  expect_equal(get_trueK(dat = epi_data2,
                         case.col = cases_n,
                         groups = c(place, country))$K,
               expt_vals2)

})
