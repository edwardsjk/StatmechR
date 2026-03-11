## MISC FXNS ##

#' Simulate artificial covariates for locatons
#'
#' @param N A numeric object providing the number of locations that should be simulated.
#' @param poverty A character object providing the level of poverty in the area. Input options are `"high"`, `"average"` or `"low"`. The default is `"average"`.
#' @param climate A numeric vector between 0 and 1 providing the proportion of locations that should be assigned to the following three climate types: "wet", "moderate" and "arid". The default is 0.33 for each climate type
#' @param temperature A numeric vector providing temperature (in celcius) associated with the three different climate types. The default is `c(30, 25, 25)`.
#' @param water_san A character object providing the level of water and sanitation in the area. Input options are `"high"`, `"average"` or `"low"`.  The default is `"average"`.
#' @param stand_water A character object providing the level of standing water in the area.
#' Input options are `"high"`, `"average"` or `"low"`.  The default is `"average"`.
#' @param education A character object providing the level of education in the area. Input options are `"high"`, `"average"` or `"low"`.  The default is `"average"`.
#'
#' @return A dataframe with a row for each location and a column for each covariate.
#'
#' @importFrom dplyr tibble
#' @importFrom dplyr mutate
#'
#' @export
#'
#' @examples
#' sim_covariates(N = 300,
#'                poverty = "high",
#'                climate = c(0.3, 0.4, 0.4),
#'                temperature = c(40, 35, 50),
#'                water_san = "low",
#'                stand_water = "average",
#'                education = "low")
#'
sim_covariates <- function(N,
                           poverty = "average",
                           climate = c(0.33, 0.33, 0.33),
                           temperature = c(30, 20, 35),
                           water_san = "average",
                           stand_water = "average",
                           education = "average") {

  watsan_loc <- ifelse(water_san == "average", 0, ifelse(water_san == "low", 1.5, -1.5))

  stand_loc <- ifelse(stand_water == "average", 0, ifelse(stand_water == "low", 1.5, -1.5))

  edu_loc <- ifelse(education == "average", 0, ifelse(education == "low", 1.5, -1.5))

  if(poverty == "average"){

    rc <- tibble(pov = rbeta(N, 3, 3, 0.5))

  }

  if(poverty == "low"){

    rc <- tibble(pov = rbeta(N, 1, 3, 0))

  }

  if(poverty == "high"){

    rc <- tibble(pov = rbeta(N, 3, 1, 0))

  }

  rc <- rc |> mutate(
    clim = sample(c("wet","mod","arid"), N, replace = T, prob = climate), #climate
    temp = rnorm(N, (clim == "wet") * temperature[1] + (clim == "mod") * temperature[2] + (clim == "arid") * temperature[3], 5), #Temprature
    watsan = plogis(rnorm(N, 2 + -2 * pov), location = watsan_loc, scale = 0.5), #wat/san access
    stand = plogis(rnorm(N, -3 + 2 * (clim == "wet")+ -2 * (clim == "dry") + 2 * (clim == "wet" | clim == "mod")), location = stand_loc, scale = 0.5),#standing water
    educ = plogis(rnorm(N, 1 + -2 * pov), location = edu_loc, scale = 0.5) #Education
  )

  return(rc)

}

#' Compute R0 values from simulated covariate data
#'
#' @param covs Dataframe providing the covariate values for each location
#' @param b0 Numeric object providing the R0 intercept on the log scale. The default is `log(2)`
#' @param sd Numeric object providing the standard deviation in R0 on the log scale. The default is `0.1`.
#' @param b1 Numeric object providing the coefficient for the wat/san access covariate. The default is `log(0.75)`.
#' @param b2 Numeric object providing the coefficient for the standing water covariate. The default is `log(1.2)`
#' @param b3 Numeric object providing the coefficient for the interaction of wat/san access and standing water covariates. The default is `log(1.2)`
#' @param temp_scale Numeric object providing a scale parameter for a logistic distribution of the temperature covariate. The default is `4`
#' @param temp_loc Numeric object providing a location parameter for a logistic distribution of the temperature covariate. The default is `20`
#' @param b4 Numeric object providing the coefficient for the temperature covariate. The default is `log(1.2)`
#'
#' @return A numeric vector providing the calculated R0 values for each location
#'
#' @importFrom dplyr mutate
#'
#' @export
#'
#' @examples
#'
#' compute_Rs(covariate_data, b0 = log(2), sd = .1,
#'            b1 = log(.75), b2 = log(1.2), b3 = log(1.2),
#'            temp_scale = 4, temp_loc = 20, b4 = log(1.2))
#'
compute_Rs <- function(covs, b0 = log(2), sd = .1, b1 = log(.75),
                       b2 = log(1.2), b3=log(1.2),
                       temp_scale = 4, temp_loc=20,
                       b4 = log(1.2) ) {

  ##first normalize the covariates
  covs <- covs |> mutate(watsan_stand = ((1-watsan)*stand-mean((1-watsan)*stand))/sd((1-watsan)*stand), #normalise interactoin term
                         watsan = (watsan-mean(watsan))/sd(watsan),
                         stand = (stand-mean(stand)/sd(stand)),
                         R0 = exp(rnorm(n(),b0+b1*watsan+b2*stand+b3*watsan_stand +
                                          b4*plogis(temp, temp_loc, temp_scale), sd)))

  return(covs$R0)

}
