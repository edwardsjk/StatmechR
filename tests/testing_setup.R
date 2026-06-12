library(tidyverse)

r_files <- list.files(path = "/users/a/b/abagaels/StatMechR/R")

r_names <- str_remove(r_files, pattern = ".R")

for(i in r_names){

  use_test(i)

}

## Test single files ##

test_file("/users/a/b/abagaels/StatMechR/tests/testthat/test-utils.R")
test_file("/users/a/b/abagaels/StatMechR/tests/testthat/test-plot_funcs.R")
test_file("/users/a/b/abagaels/StatMechR/tests/testthat/test-error_funcs.R")



