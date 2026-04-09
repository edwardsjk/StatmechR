#### PACKAGE DEVELOPMENT: STATMECHR ####

library(roxygen2)
library(testthat)
library(knitr)
library(checkmate)
library(cli)

#create_package("/users/a/b/abagaels/StatMechR")

## ??

devtools::dev_sitrep()
devtools::update_packages("devtools")
devtools::install_dev_deps()

setwd("/users/a/b/abagaels/StatMechR")
document()

## Add build ignore file

usethis::use_build_ignore(c("build", "data-raw", "old_funcs", ".Rprofile"), escape = T)

## Add dependencies

packages <- c("dplyr", "lubridate", "ggplot2", "zoo", "doParallel", "foreach",
              "SuperLearner", "gam", "rpart", "randomForest", "e1071", "stringr",
              "parallel")

for(i in packages){

  use_package(i, type = "Imports")

}

document()

## Add vignettes

usethis::use_vignette("Vignette_1")


## Check package

devtools::check()



