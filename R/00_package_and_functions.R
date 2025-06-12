
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
#' @param ci_type character - The type of confidence interval calculatgion for P and Q matrices can be either normal, none, or bootstrap
#' @param n_cores numeric Indicates how many cores we want to use for the bootstrapping ci. Currently only works with 1 core
#' @param B numeric Indicates number of bootstraps
#' @param conf_level numeric - confidence intervals level to be used
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
                           calc_spline   = TRUE,
                           ci_type       = "bootstrap",
                           n_cores       = 10,
                           B             = 2,
                           conf_level    = 0.95,
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
    spline_basis <-  predict(spline_basis_fit,
                             newx = prediction_grid$age) |>
      as_tibble() |>
      set_names(paste0("age_spline", 1:spline_df))
    
    prediction_grid <- bind_cols(prediction_grid, spline_basis)
    
    spline_basis_fit <- spline_basis_fit |>
      as_tibble() |>
      set_names(paste0("age_spline", 1:spline_df))
    
    result <- .data |>
      bind_cols(spline_basis_fit)
    
  } else {
    
    spline_basis <- prediction_grid |> 
      dplyr::select(age)
    
    result <- .data
    
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
  
  # prepare spline basis
  
  # calculate model
  result <- result |>
    group_by(across(all_of(strat_vars))) |>
    group_nest() |>
    # just in case
    na.omit() |>
    # removing solitary observations
    mutate(data = map(data, ~ .x |> 
                        group_by(hhidpn) |> 
                        # supposed to remove solitary observations.
                        filter(n() > 1) |> 
                        ungroup())) |>
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
  # calculate the rate and prob. lists with CI
  # if none - simply calculate the matrices and return them
  if(ci_type == "none") { 
    
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
    
    # add ci and done
    # this is the final result
    result_df <- rate |>
      # join together
      full_join(prob) |>
      # remove D-D transitions (empty)
      filter(from != "State 3") |>
      # remove recovery possibility (empty)
      filter(!(from == "State 2" & to == "State 1"))
    
  } else { 
    
    # make data with explicit variables for list calculation
    pred <- prediction_grid |>
      # add our age interval
      mutate(age_interval = age_int) |>
      # join the model
      left_join(result, by = strat_vars) |>
      mutate(ci    = ci_type,
             cl    = conf_level,
             B     = B,
             cores = n_cores) |>
      mutate(Q = list(Q)) 
    # |>
    #   # here for test and speed we can take a small portion of data
    #   slice(1:2)  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # 
    # rate list first
    pred$q_list <- pmap(
      list(pred$model,
           pred$data_fit,
           pred$ci,
           pred$cl,
           pred$B,
           pred$cores),
      ~ qmatrix.msm(x          = ..1, 
                    covariates = ..2,
                    ci         = ..3,
                    cl         = ..4,
                    B          = ..5,
                    cores      = ..6)
    )
    
    # probability second
    pred$p_list <- pmap(
      list(pred$model,
           pred$age_interval,
           pred$data_fit,
           pred$ci,
           pred$cl,
           pred$B,
           pred$cores),
      ~ pmatrix.msm(x          = ..1, 
                    t          = ..2, 
                    covariates = ..3,
                    ci         = ..4,
                    cl         = ..5,
                    B          = ..6,
                    cores      = ..7)
    )
    
    # --------------------------------------------------------------#
    # PT4 tidy the reslting lists and return a nice dataframe in the end
    # from here to the end only data manipulation to retrieve the final 
    # dataset with columns lower, upper and estimate
    # nothing new is calculated
    # tidy the q_list
    
    result_df <- pred |>
      mutate(
        q_tidy = map(q_list, ~ {
          # extract components
          lower <- .x$L |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from,
                         names_to  = "to",
                         values_to = "lower")
          
          estimate  <- .x$estimate |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from,
                         names_to  = "to",
                         values_to = "estimate")
          
          upper <- .x$U |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from,
                         names_to  = "to",
                         values_to = "upper")
          
          lower |>
            left_join(estimate, by = c("from", "to")) |>
            left_join(upper, by = c("from", "to")) |>
            mutate(type = "q")
        }
        )) |>
      # now we do the same routine for tidying p
      mutate(
        p_tidy = map(p_list, ~ {
          # extract components
          lower <- .x$L |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from,
                         names_to  = "to",
                         values_to = "lower") |>
            mutate(from = str_c("State", from, sep = " "),
                   to   = str_c("State", parse_number(to), sep = " "))
          
          estimate  <- .x$estimates |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from,
                         names_to  = "to",
                         values_to = "estimate")
          
          upper <- .x$U |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from,
                         names_to  = "to",
                         values_to = "upper") |>
            mutate(from = str_c("State", from, sep = " "),
                   to   = str_c("State", parse_number(to), sep = " "))
          
          lower |>
            left_join(estimate, by = c("from", "to")) |>
            left_join(upper, by = c("from", "to")) |>
            mutate(type = "p")
        }
        )) |>
      # remove helper columns combine P and Q tidy data and unnest them
      # The column type now indicates p or q
      mutate(
        combined_tidy = map2(q_tidy, p_tidy, bind_rows)) |>
      # there are helper columns that we do not need no more
      select(-c(data_fit, age_interval, data, model,
                q_list, p_list, q_tidy, p_tidy, ci, cl, B, cores, Q)) |>  
      unnest(combined_tidy) |>
      # remove D-D transitions (empty)
      filter(from != "State 3") |>
      # remove recovery possibility (empty)
      filter(!(from == "State 2" & to == "State 1"))
    
  }

  # return the final dataframe
  
  return(result_df)
  
}
# ------------------------------------------------------------------- #