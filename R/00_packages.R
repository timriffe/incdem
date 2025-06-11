
# whenever we add packages, after installing locally, run renv::snapshot()
# renv::snapshot()

# renv::restore()
# ------------------------------------------------------------------- #
library(tidyverse)
library(splines)
library(msm)
library(haven)
library(zoo)
library(slider)
library(parallel)
library(doParallel)
# ------parallel# ------------------------------------------------------------------- #
# impute age using some simple logic
impute_age <- function(age, wave){
  if (all(is.na(age))){
    return(age)
  }
  age <-  zoo::na.approx(age, x = wave, na.rm = FALSE)
  if (any(is.na(age))){
    if (is.na(age[1])){
      # handle leading NAs
      diffs            <- diff(wave) * 2
      nNAs             <- rle(is.na(age))$lengths[1]
      first_non_NA_age <- age[!is.na(age)][1]
      subtract_this    <- diffs[1:nNAs] |> rev() |> cumsum() |> rev()
      age[1:nNAs]      <- first_non_NA_age - subtract_this
    }
    n <- length(age)
    if (is.na(age[n])){
      diffs            <- diff(wave) * 2
      nNAs             <- rle(is.na(rev(age)))$lengths[1]
      last_non_NA_age  <- age[!is.na(age)] |> rev() %>% '['(1)
      add_this         <- diffs[(n-nNAs):(n-1)] |> cumsum()
      age[(n-nNAs+1):n]  <- last_non_NA_age + add_this
    }
  }
  age
}


#' Fit Multistate Model with Sensitivity to Covariates
#'
#' Fits a multistate model by including one or more covariates and generating predicted transition rates and probabilities.
#' The function nests the data by stratification and covariate levels, fits the MSM model for each group, 
#' then produces a prediction grid for each combination of covariates and age, and computes transition rates and probabilities.
#'
#' @param .data A dataframe with individual data. in our case hrs_msm
#' @param strat_vars A character vector of variable names to stratify the MSM (group and nest)
#' @param covariate_var A character vector of covariate variable names to include in the model.
#' @param age_int A numeric value specifying the interval for prediction ages (default is 0.25 years).
#' @param cont_grid This parameter is only used if one of the stratification or vovariate vars is continuous. Say we have continuous time, in this case this variable should take the values like 2000, 2010, for data fitting
#' @return A data frame with estimated transition rates and probabilities 
#' 

# ------------------------------------------------------------------- #
# check the comments and example code below
# this function has 4 (four) parts:
# 1 creating the grid, 
# 2 fitting the msm
# 3 predicting from the created grid and model
# 4 reshaping and returning the outcome
# can take any variables as covariate_var or strat_vars
# if covariate or strat var is continuous cont_grid should be provided
# it should be something ordinal or nominal like c(2000, 2010) 
# a assume that age is always provided 
# and 3 splines exist as covariates by default for now 
# NOTE: Dont pay attention to warning. IT DOES NOT AFFECT RESULTS
# I suggest using newly created obs_date as a time variable
fit_msm_sensitivity <- function(.data,
                                strat_vars    = NULL,
                                covariate_var = NULL,
                                age_int       = 0.25,
                                cont_grid     = NULL) {
  # --------------------------------------------------------------#
  # PT1
  # We start by constructing the grid
  # the ifelse part takes arbitrary specified columns and 
  # makes a list with the correspodning unique values
  # This first part is only used in case we have continuous covariate 
  if(!is.null(cont_grid)) {
    
    vars        <- map(c(strat_vars), 
                       ~ unique(na.omit(.data[[.x]])))
    
    vars[[(length(vars) + 1)]] <- cont_grid 
    
    names(vars) <- c(strat_vars, covariate_var)
    
    # othervise we use this part
  } else { 
    
    vars        <- map(c(covariate_var, strat_vars), 
                       ~ unique(na.omit(.data[[.x]]))) |>
      set_names(c(covariate_var, strat_vars))
    
  }
  
  # add base predictors. in our model it is age
  base_grid <- list(
    age    = seq(50, 100, by = age_int)
  )
  
  # construct prediction grid dynamically
  prediction_grid <- exec(crossing, !!!c(base_grid, vars))
  
  # predict the spline basis for the age specification
  spline_basis <- predict(spline_basis_fit,
                          newx = prediction_grid$age) |>
    as.data.frame() |>
    setNames(paste0("age_spline", 1:3))
  
  # bind data together
  prediction_grid <- bind_cols(prediction_grid, spline_basis) |>
    # mutate(female = ifelse(female == 1, "Female", "Male")) |>
    # nest data by correspodnding covariates
    group_by(across(all_of(c(
      strat_vars, covariate_var, "age"
    )))) |>
    group_nest(keep = TRUE, .key = "data_fit") |>
    # create a nested list of fitting values
    mutate(data_fit = map(data_fit, ~ .x |>
                            dplyr::select(-c(
                              all_of(strat_vars), age
                            )) |>
                            as.list()))
  
  # --------------------------------------------------------------#
  # PT2
  # create fit model with arbitrary covariates
  # always assume that age aplines are used and age is provided
  # NOTE: Takes time
  base_covs   <- c("age_spline1", "age_spline2", "age_spline3")
  all_covs    <- c(base_covs, covariate_var)
  cov_formula <- reformulate(all_covs)
  
  # Stratify by one or more variables
  result <- .data |>
    group_by(across(all_of(strat_vars))) |>
    group_nest() |>
    na.omit() |>
    mutate(model = map(
      data,
      ~ msm(
        formula    = state_msm ~ age,
        subject    = hhidpn,
        data       = .x,
        qmatrix    = Q,
        obstype    = obstype,
        deathexact = 3,
        covariates = cov_formula,
        control    = list(fnscale = 5000, maxit = 25000),
        gen.inits  = TRUE,
        method     = "BFGS"
      )
    )) |>
    dplyr::select(-data)
  
  # --------------------------------------------------------------#
  # PT3
  # here we predict the model by the corresponding created grid and 
  # calculate the rate and prob. matrices
  pred <- prediction_grid |>
    # add our age interval
    mutate(age_interval = age_int) |>
    # join the model
    left_join(result, by = strat_vars) |>
    # calculate rate
    mutate(
      rate = map2(
        .x = model,
        .y = data_fit,
        ~ qmatrix.msm(x          = .x, 
                      covariates = .y)$estimates |>
          as.table() |>
          as.data.frame() |>
          rename(from = Var1, 
                 to   = Var2, 
                 rate = Freq)
      ),
      # calculate probability
      prob = pmap(
        list(model, data_fit, age_interval),
        ~ pmatrix.msm(x          = ..1, 
                      t          = ..3, 
                      covariates = ..2) |>
          as.table() |>
          as.data.frame() |>
          rename(from = Var1, 
                 to   = Var2, 
                 prob = Freq)
      )
    )
  # ------------------------------------------------------------------- #
  # PT4 
  # return resulting data with rates and probabilities 
  # unnest rate and chosen variables
  rate <- pred |>
    dplyr::select(all_of(c(strat_vars, covariate_var)), age, rate) |>
    unnest(rate)
  
  # unnest probability
  prob <- pred |>
    dplyr::select(all_of(c(strat_vars, covariate_var)), age, prob) |>
    unnest(prob)
  
  # this is the final result
  result_df <- rate |>
    # join together
    full_join(prob) |>
    # remove D-D transitions (empty)
    filter(from != "State 3") |>
    # remove recovery possibility (empty)
    filter(!(from == "State 2" & to == "State 1"))
  
  return(result_df)
  
}