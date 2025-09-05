library("devtools")
library(roxygen2)
setwd("/users/a/b/abagaels/StatMechR")
#document()

## Add dependencies

packages <- c("dplyr", "lubridate", "ggplot2", "zoo", "doParallel", "foreach",
              "SuperLearner", "gam", "rpart", "randomForest", "e1071",
              "bartMachine", "stringr", "parallel")

for(i in packages){

  use_package(i, type = "Imports")

}

document()







