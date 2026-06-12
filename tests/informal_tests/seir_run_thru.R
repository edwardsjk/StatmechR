#### CHOLERA SEIR ####

devtools::load_all()

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

alldat <- read.csv("/users/a/b/abagaels/Data/statmech_data/yemen_cholera/df_cases_updated.csv")

alldat$n_cases <- as.numeric(ifelse(is.na(alldat$n_cases), 0, alldat$n_cases))
alldat$epiweek_date <- as.Date(alldat$epiweek_date)

alldat <- alldat %>%
  dplyr::filter(epiweek_date <= as.Date("2018-01-01"))

alldat$district[which(alldat$governorate == "al_bayda" & alldat$district == "az_zahir")] <- "az_zahir(ab)"
alldat$district[which(alldat$governorate == "al_jawf" & alldat$district == "az_zahir")] <- "az_zahir(aj)"

tmp <- alldat %>%
  group_by(governorate, district) %>%
  dplyr::slice(1) %>%  #filtered to just the first infection in each district
  arrange(district) %>%
  ungroup()

population <- tmp$pop_2017

trueK <- get_trueK(alldat, case.col = n_cases, groups = c(governorate, district)) %>%
  arrange(district)

plot_true_curves(alldat, X = epiweek_date, plot.group = district, count = n_cases, legend = F)

mask_dat <- alldat %>%
  filter((as.Date(epiweek_date) <= as.Date("2017-07-16")) & is.na(epiweek_date) == F) %>%
  group_by(district) %>%
  arrange(epiweek_date) %>%
  mutate(t = 1:n(), first_week = which(n_cases > 0)[1],
         epi_t = t - (first_week - 1), c_cases = cumsum(n_cases),
         pop_dens = pop_2017/area_sqkm) %>%
  ungroup() %>%
  arrange(district, epiweek_date)

plot_mask_curves(dat = alldat, maskdat = mask_dat, X = epiweek_date, count = n_cases,  plot.group = governorate)

mask_dat <- mask_dat %>%
  select(-governorate, -epiweek_date, -t, -first_week, -code, -area_sqkm) %>%
  mutate(epi_t = ifelse(is.na(epi_t) == T, 0, epi_t))

mask_dat$pop_type[which(mask_dat$pop_type == "Rural")] <- 1
mask_dat$pop_type[which(mask_dat$pop_type == "Difficult")] <- 2
mask_dat$pop_type[which(mask_dat$pop_type == "Very difficult")] <- 3
mask_dat$pop_type[which(mask_dat$pop_type == "Urban")] <- 4

mask_dat$camp[which(mask_dat$camp == "houti")] <- 1
mask_dat$camp[which(mask_dat$camp == "coalition")] <- 2

mask_dat[, c(2:54)] <- as.data.frame(sapply(mask_dat[, c(2:54)], as.numeric))

na_cols <- names(which(sapply(mask_dat, anyNA)))

ecs <- split(mask_dat$n_cases, f = mask_dat$district)

ecs2 <- lapply(ecs, function(x){

  if(any(x > 0) == TRUE){

    new_ec <- x[which(x > 0)[1]:length(x)]

    if(length(new_ec) > 3){

      sm_line <- round(smooth.spline(x = 1:length(new_ec), y = new_ec, spar = 0.5)$y)
      sm_line[which(sm_line < 0)] <- 0

      sm_line

    }else{

      new_ec

    }

  }else{

    0

  }

})

enddate <- as.Date(max(alldat$epiweek_date, na.rm = T))
tau <- round(as.numeric(as.Date(max(alldat$epiweek_date, na.rm = T)) - as.Date("2017-07-16"))/7)

curcases <- sapply(ecs, sum)

#### Epi alone ####

# starting_vals <- data.frame(I0 = rep(5, length(ecs)),
#                             beta = rep(0.1, length(ecs)),
#                             prob_detect = rep(0.5, length(ecs)))

starting_vals <- data.frame(I0 = rep(5, length(ecs)),
                            beta = rep(0.1, length(ecs)))

epimdl <- run_custom_model(ecs = ecs2, N = population,
                           epi.mdl.func = yemen_model_new,
                           epi.mdl.pars = starting_vals,
                           error.func = poisson_error,
                           tau = tau, timestep = 7)

k_mech <- as.vector(do.call(rbind, lapply(epimdl, function(x){x$estK})))

model_params <- as.data.frame(do.call(rbind, lapply(epimdl, function(x){x$params})))

cov_dat <- mask_dat %>%
  distinct(.keep_all = T) %>%
  filter(epi_t >= 0) %>%
  mutate(c_cases_1000 = c_cases/(pop_2017/1000))

xvars_wt <- cov_dat %>%
  select(-all_of(na_cols), -district, -n_cases, -c_cases, -c_cases_1000, -pop_2017)

outcome = cov_dat$c_cases_1000

mod <- SL_fit(x = xvars_wt, y = outcome, family = "gaussian")

xvars_pred <- mask_dat %>%
  distinct(.keep_all = T) %>%
  group_by(district) %>%
  filter(epi_t == max(epi_t)) %>%
  ungroup() %>%
  select(-district, -n_cases, -c_cases, -pop_2017) %>%
  mutate(epi_t = epi_t + tau) %>%
  select(-all_of(na_cols))

k_stat <- SL_pred(mdl = mod, x = xvars_pred)

k_stat2 <- k_stat * (population/1000)

comb_vars <- cov_dat %>%
  group_by(district) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  select(-all_of(na_cols), -n_cases, -c_cases, -epi_t, -district,
         -pop_2017, -c_cases_1000)

# starting_vals <- data.frame(I0 = rep(5, length(ecs)),
#                             beta = rep(0.1, length(ecs)),
#                             prob_detect = rep(0.5, length(ecs)))

starting_vals <- data.frame(I0 = rep(5, length(ecs)),
                            beta = rep(0.1, length(ecs)))

k_combinf <- run_combined_model(epi.curves = ecs2, covdat = comb_vars,
                                pop.N = population, stat.model = "SL",
                                epi.model = "custom",
                                custom.model = yemen_model_new,
                                starting.vals = starting_vals,
                                error.func = poisson_error,
                                penalty.func = poispen, tau = tau,
                                timestep = 7, max.iter = 5,
                                epi.parallel = F, cores = NULL)

k_iter <- k_combinf$Khist

bias_plot <- plot_bias(trueK, k_iter)
