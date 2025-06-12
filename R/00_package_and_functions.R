
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
library(furrr)
library(future)
library(multidplyr)
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
#' @param age_from_to A numberic vector. This parameter controls the lower and upper bnound of the age prediction grid 
#' @param spline_df A number. This parameter indicates the degrees of freedom for the age spline fitting. E.G. if 3 then 3 splines are created 
#' @param spline_type A character. Either "ns" or "bs" - indicates the type of spline to be fitted to the age
#' @param calc_spline Logical. Indicates weather we want to calculate spline TRUE or not FALSE
#' @param n_cores numeric Indicates how many cores we want to use for the bootstrapping ci. Currently only works with 1 core
#' @param B numeric Indicates number of bootstraps
#' @param Q matrix Indicates a Q mateix that will be used for model fitting. Diagonal elemens well be calculated automatically

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
fit_msm_models <- function(.data,
                                # main variables
                                strat_vars    = NULL,
                                covariate_var = NULL,
                                # only used for continuous variables
                                cont_grid     = NULL,
                                # age grid specification part
                                age_from_to   = c(50, 100),
                                age_int       = 0.25,
                                # spline specification part
                                spline_df     = NULL,
                                # do we want to calculate spline at all?
                                spline_type   = "ns",
                                calc_spline   = FALSE,
                                n_cores       = 1,
                                B             = 2,
                                # create Q matrix
                                Q = rbind(
                                  c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
                                  c(0, -.01,   0.1),  # dementia can go to death
                                  c(0, 0,   0)        # death is absorbing
                                )) {

  
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
  # from age_pred_grid from fisr to last by age_int
  base_grid <- list(
    age    = seq(age_from_to[1], 
                 age_from_to[2], 
                 by = age_int)
  )
  
  # construct prediction grid dynamically
  prediction_grid <- exec(crossing, !!!c(base_grid, vars))
  
  if(calc_spline) { 
    
    # construct spline basis for age
    # Here we can also add another types of splines
    # currently only ns or bs
    # spline_df specified number of splines
    spline_basis_fit <- get(spline_type)(.data$age, df = spline_df)
    
    # predict the spline basis for the age specification
    spline_basis <- predict(spline_basis_fit,
                            newx = prediction_grid$age) |>
      as_tibble() |>
      set_names(paste0("age_spline", 1:spline_df))
    
    prediction_grid <- bind_cols(prediction_grid, spline_basis)
    
  } else { 
    
    spline_basis <- prediction_grid |> 
      dplyr::select(age)
    
  }
  
  # bind data together
  prediction_grid <- prediction_grid |>
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
  # always assume that age splines are used and age is provided
  # NOTE: Takes time
  base_covs   <- names(spline_basis)
  all_covs    <- c(base_covs, covariate_var)
  cov_formula <- reformulate(all_covs)
  
  # calculate model
  result <- .data |>
    group_by(across(all_of(strat_vars))) |>
    group_nest() |>
    na.omit() |>
    mutate(model = map(data,
                       ~ msm(
                         formula    = state_msm ~ age,
                         subject    = hhidpn,
                         data       = .x,
                         qmatrix    = Q,
                         obstype    = obstype,
                         deathexact = 3,
                         covariates = cov_formula,
                         control    = list(fnscale = 5000, 
                                           maxit   = 25000),
                         gen.inits  = TRUE,
                         method     = "BFGS"
                       )
    ))
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
  
  # --------------------------------------------------------------#
  # PART 4
  # boothstrapping for confidence intervals
  # ci <- prediction_grid |>
  #   mutate(age_interval = age_int) |>
  #   left_join(result, by = strat_vars) |>
  #   # slice(1:2) |>
  #   mutate(
  #     q_list = map2(model,
  #                   data_fit, ~
  #                     boot_qmatrix(.x, .y,
  #                                  B     = B,
  #                                  # not working with > 1
  #                                  cores = n_cores))
  #   ) |>
  #   # now calculate the confidence intervals from the bootstraps
  #   # and unnest
  #   mutate(
  #     q_array = map(q_list, simplify2array),
  #     q_ci = map(q_array, ~ list(
  #       lower = apply(.x, c(1, 2), quantile, probs = 0.025, na.rm = TRUE),
  #       upper = apply(.x, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
  #     )),
  #     q_tidy = map(q_ci, ~ {
  #       # extract components
  #       lower_mat <- .x$lower
  #       upper_mat <- .x$upper
  #       
  #       # convert each to data frames
  #       lower_df <- as.data.frame(as.table(lower_mat)) |>
  #         set_names(c("from", "to", "lower"))
  #       upper_df <- as.data.frame(as.table(upper_mat)) |>
  #         set_names(c("from", "to", "upper"))
  #       
  #       lower_df |>
  #         left_join(upper_df, by = c("from", "to")) |>
  #         mutate(type = "q")
  #     }
  #     ))|>
  #   dplyr::select(-c(data_fit, data, model, 
  #                    age_interval, q_array, q_ci, rate)) |>
  #   # same can be done for p matrix
  #   # Select only what you need and unnest
  #   unnest(q_tidy)
  
  
  
  # ------------------------------------------------------------------- #
  # PT5
  # return resulting data with rates and probabilities 
  # unnest rate and chosen variables
  rate <- pred |>
    dplyr::select(all_of(c(strat_vars, covariate_var)), age, rate) |>
    unnest(rate)
  
  # unnest probability
  prob <- pred |>
    dplyr::select(all_of(c(strat_vars, covariate_var)), age, prob) |>
    unnest(prob)
  
  # add ci and done
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

# ------------------------------------------------------------------- #
boot_qmatrix <- function(fitted_model,
                         new_covariates,
                         B     = 1,
                         cores = 1) {
  stat_fn <- function(m) qmatrix.msm(m, covariates = new_covariates)$estimates
  boot.msm(fitted_model,
           stat  = stat_fn,
           B     = B,
           cores = cores)
}
