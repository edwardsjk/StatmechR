#### PACKAGE DEVELOPMENT: STATMECHR ####

library(devtools)
library(roxygen2)
library(testthat)
library(knitr)

#create_package("/users/a/b/abagaels/StatMechR")

devtools::dev_sitrep()
devtools::update_packages("devtools")
devtools::install_dev_deps()

setwd("/users/a/b/abagaels/StatMechR")
document()

usethis::use_build_ignore(c("build", "data-raw", "funcs_in_progress"), escape = T)

## Add dependencies

# packages <- c("dplyr", "lubridate", "ggplot2", "zoo", "doParallel", "foreach",
#               "SuperLearner", "gam", "rpart", "randomForest", "e1071",
#               "bartMachine", "stringr", "parallel")
#
# for(i in packages){
#
#   use_package(i, type = "Imports")
#
# }
#
# document()

## Add vignettes

usethis::use_vignette("Vignette_1")





