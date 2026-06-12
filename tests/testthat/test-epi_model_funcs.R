## TEST EPI_MODEL_FUNCS.R ##

devtools::load_all()

########################
#### gaussian_mod() ####
########################

test_that("gaussian_mod() gives an error if there are not enough parameter values", {

  pars_1 <- c(1, 3)

  expect_error(gaussian_mod(pars = pars_1, times = 10),
               "gaussian_mod\\(\\) requires 3 parameter values")

})

test_that("gaussian_mod() gives an error if there are not enough parameter values", {

  pars_1 <- c(1, 3, 0.1)

  expect_equal(length(gaussian_mod(pars = pars_1, times = 10)),
               10)

})


##############################
#### run_gaussian_model() ####
##############################

params <- data.frame(estK = sample(50:200, 20, replace = T),
                     peak = sample(5:15, 20, replace = T),
                     spread = sample(seq(0.5, 1.5, 0.05), 20, replace = T))

cur_curves <- lapply(1:20, function(x){

  round(gaussian_mod(pars = params[x,], times = 20))

})

starting_vals <- data.frame(estK = rep(100, times = 20),
                            peak = rep(10, times = 20),
                            spread = rep(1, times = 20))

test_that("run_gaussian_model() returns object of correct length", {

  tmp <- run_gaussian_model(ecs = cur_curves, epi.mdl.pars = starting_vals,
                            error.func = poisson_error)

  expect_equal(length(tmp), length(cur_curves))

})

test_that("run_gaussian_model() fitted curve closely matches true curves", {

  tmp <- run_gaussian_model(ecs = cur_curves, epi.mdl.pars = starting_vals,
                            error.func = poisson_error)

  for(i in length(tmp)){

    true_curve <- cur_curves[[i]]

    fitted_pars <- tmp[[i]][[1]]

    fitted_curve <- gaussian_mod(pars = fitted_pars, times = 20)

    expect_gt(cor(fitted_curve, true_curve), 0.99)

    model_rmse <- sqrt(mean((fitted_curve - true_curve)^2))
    naive_rmse <- sqrt(mean((mean(true_curve) - true_curve)^2))
    expect_lt(model_rmse, naive_rmse * 0.05)

  }

})

############################
#### run_custom_model() ####
############################

yemen_model_new <- function(N, pars, time_step, tau){

  #pars includes:
  #I0, beta, prob_detect

  phi = 0.2
  #prob_detect = exp(-abs(log(abs(pars[3]))))
  #prob_detect = abs(pars[3])
  prob_detect = 0.425

  I0 = round(pmin(abs(pars[1]), N))

  if(I0 < 2){

    I0 = 2

  }

  beta = abs(pars[2])

  cur_state <- c(t = 0, S = N - I0,
                 I = I0, R = 0,
                 incident = round(I0 * prob_detect))

  tmp <- list(cur_state)
  i <- 1

  while(i <= round(tau)){

    i = i + 1

    lambda <- beta * (cur_state[3]/N)

    StoI = cur_state[2] * (1-exp(-time_step*lambda))
    StoI = ifelse(is.na(StoI), 0, round(StoI))

    ItoR = cur_state[3] * (1-exp(-time_step*phi))
    ItoR = ifelse(is.na(ItoR), 0, round(ItoR))

    cur_state[1] <- cur_state[1] + time_step
    cur_state[2] <- cur_state[2] - StoI
    cur_state[3] <- cur_state[3] + StoI - ItoR
    cur_state[4] <- cur_state[4] + ItoR
    cur_state[5] <- round(StoI * prob_detect)

    tmp[[i]] <- cur_state

  }

  epi <- bind_rows(tmp)

  return(epi)

}

# params2 <- data.frame(I0 = sample(1:10, 20, replace = T),
#                       beta = sample(seq(0.1, 0.2, 0.01), 20, replace = T),
#                       prob_detect = sample(seq(0.25, 0.75, 0.05), 20, replace = T))

params2 <- data.frame(I0 = sample(1:10, 20, replace = T),
                      beta = sample(seq(0.1, 0.2, 0.01), 20, replace = T))

population <- sample(1000:10000, 20)

cur_curves2 <- lapply(1:20, function(x){

  yemen_model_new(population[x], unname(unlist(params2[x,])), time_step = 7, tau = 12)$incident

})

# starting_vals2 <- data.frame(I0 = rep(5, times = 20),
#                              beta = rep(0.15, times = 20),
#                              prob_detect = rep(0.5, times = 20))

starting_vals2 <- data.frame(I0 = rep(5, times = 20),
                             beta = rep(0.15, times = 20))

test_that("run_custom_model() closely matches the fitted curve",{

  tmp <- run_custom_model(ecs = cur_curves2, N = population,
                          epi.mdl.func = yemen_model_new,
                          epi.mdl.pars = starting_vals2,
                          error.func = poisson_error2,
                          tau = 12,
                          timestep = 7)

  for(i in length(tmp)){

    true_curve2 <- cur_curves2[[i]]

    fitted_curve2 <- tmp[[i]][[4]][1:length(true_curve2)]

    expect_gt(cor(fitted_curve2, true_curve2, use = "complete.obs"), 0.99)

    model_rmse2 <- sqrt(mean((fitted_curve2 - true_curve2)^2))
    naive_rmse2 <- sqrt(mean((mean(true_curve2) - true_curve2)^2))
    expect_lt(model_rmse2, naive_rmse2 * 0.05)

  }

})

#####################
#### build_SIR() ####
#####################

######################
#### build_SEIR() ####
######################
