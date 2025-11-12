
# whenever we add packages, after installing locally, run renv::snapshot()
# renv::snapshot(force=TRUE)


# renv::restore()
# ------------------------------------------------------------------- #
library(doParallel)
library(tidyverse)
library(parallel)
library(splines)
library(slider)
library(haven)
library(expm)
library(msm)
library(zoo)
library(rsample)
# these 4 can be commented out !!!!!!
library(multidplyr)
library(future)
library(furrr)
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

# this function possibly deprecate, but it's a nice normal solution.
qmatrix.msm_wrapper <- function(x, ci = c("none", "normal", "delta", "bootstrap"),
                           B = 1000, cores = 1, age_interval = .25, ...) {
  ci <- match.arg(ci)
  Q <- NULL
  
  if (ci == "normal" && B > 1) {
    coef_hat <- coef(x)
    vcov_hat <- vcov(x)
    sim_pars <- mvtnorm::rmvnorm(B, mean = coef_hat, sigma = vcov_hat)
    
    qmat_from_par <- function(pvec) {
      x_tmp <- x
      x_tmp$paramdata$opt$par <- pvec
      suppressWarnings(qmatrix.msm(x_tmp, ci = "none", ...))
    }
    
    if (cores > 1) {
      if (.Platform$OS.type == "unix") {
        # this runs more efficiently on linux machines... like mine!
        qlist <- parallel::mclapply(1:B, function(i) qmat_from_par(sim_pars[i, ], x, ...),
                                    mc.cores = cores)
      } else {
        # this is for Windows machines
        cl <- parallel::makeCluster(cores)
        doParallel::registerDoParallel(cl)
        qlist <- foreach::foreach(i = 1:B, .packages = "msm") %dopar% {
          qmat_from_par(sim_pars[i, ], x)
        }
        parallel::stopCluster(cl)
      }
    } else {
      qlist <- lapply(1:B, function(i) qmat_from_par(sim_pars[i, ]))
    }
    
    qarray <- simplify2array(qlist)
    Q <- list(
      estimate = apply(qarray, 1:2, mean, na.rm = TRUE),
      L = apply(qarray, 1:2, quantile, probs = 0.025, na.rm = TRUE),
      U = apply(qarray, 1:2, quantile, probs = 0.975, na.rm = TRUE)
    )
    
  } else {
    Q_res <- qmatrix.msm(x, ci = ci, B = B, cores = cores, ...)
    if (ci %in% c("delta", "bootstrap", "normal")) {
      Q <- list(
        estimate = Q_res$estimates,
        L = Q_res$L,
        U = Q_res$U
      )
    } else {
      Q <- list(estimate = Q_res)
    }
  }
  
  # Always return P also
  if (!is.null(Q$L) && !is.null(Q$U)) {
    P <- list(
      estimate = expm::expm(Q$estimate * age_interval),
      L = expm::expm(Q$L * age_interval),
      U = expm::expm(Q$U * age_interval)
    )
  } else {
    P <- list(estimate = expm::expm(Q$estimate * age_interval))
  }
  
  list(Q = Q, P = P)
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
#' @param cont_grid This parameter is only used if one of the stratification or covariate vars is continuous. Say we have continuous time, in this case this variable should take the values like 2000, 2010, for data fitting
#' @param age_from_to A numeric vector. This parameter controls the lower and upper bound of the age prediction grid 
#' @param spline_df A number. This parameter indicates the degrees of freedom for the age spline fitting. E.G. if 3 then 3 splines are created 
#' @param spline_type A character. Either "ns" or "bs" - indicates the type of spline to be fitted to the age
#' @param calc_spline Logical. Indicates weather we want to calculate spline TRUE or not FALSE
#' @param ci_type character - The type of confidence interval calculation for P and Q matrices can be either normal, none, or bootstrap
#' @param n_cores numeric Indicates how many cores we want to use for the bootstrapping ci. Currently only works with 1 core
#' @param B numeric Indicates number of bootstraps
#' @param conf_level numeric - confidence intervals level to be used
#' @param Q matrix Indicates a Q matrix that will be used for model fitting. Diagonal elements well be calculated automatically

#' @return A data frame with estimated transition rates and probabilities 
#' 
# ------------------------------------------------------------------- #
# check the comments and example code below
# this function has 4 (four) parts:
# 1 Creating the FIT and PREDICT DATASETS
# 2 Estimating the MSM model and prepare the basis for prediction data
# 3 assign variables to global env
# 4 Ccalculate corresponding CI
# can take any variables as covariate_var or strat_vars
# if covariate or strat var is continuous cont_grid should be provided
# it should be something ordinal or nominal like c(2000, 2010) 
# I assume that age variable is always available in the dataset 
# I suggest using newly created obs_date as a time variable
 fit_msm <- function(.data, # .data for piping
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
   
   # ------------------------------------------------------------------- #
   # PT1 - Creating the FIT and PREDICT DATASETS
   # Allows splines and no splines.
   # Allows the use of continuous covariate like year e.g. (2010, 2012, 2020)
   if(!is.null(cont_grid)) {
     
     vars <- map(c(strat_vars), ~ unique(na.omit(.data[[.x]])))
     vars[[(length(vars) + 1)]] <- cont_grid
     names(vars) <- c(strat_vars, covariate_var)
     
     # otherwise we use this part
   } else { 
     
     vars <- map(c(covariate_var, strat_vars),~ unique(na.omit(.data[[.x]]))) |>
       set_names(c(covariate_var, strat_vars))
     
   }
   
   # add base predictors. in our model it is age
   # from age_pred_grid from fisr to last by age_int
   base_grid <- list(
     age = seq(age_from_to[1], age_from_to[2], by = age_int)
   )
   
   # construct prediction grid dynamically
   prediction_grid <- exec(crossing, !!!c(base_grid, vars))
   
   if(calc_spline) { 
     
     # construct spline basis for age
     # Here we can also add another types of splines in future
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
     
     # bind data together
     prediction_grid <- prediction_grid |>
       # nest data by corresponding covariates
       group_by(across(all_of(c(strat_vars, covariate_var, "age")))) |>
       group_nest(keep = TRUE, .key = "data_fit") |>
       # create a nested list of fitting values
       mutate(data_fit = map(data_fit, ~ .x |>
                               dplyr::select(-c(all_of(strat_vars), age)) |>
                               as.list())
              )
     
   } else { # if we do not want to use splines
     # The difference here is in the data_fit. Previous removes age, this keeps it
     
     spline_basis <- prediction_grid |>
       dplyr::select(age)
     
     # here we assign original data for compatibility with other method
     result <- .data
     
     # bind data together
     prediction_grid <- prediction_grid |>
       # nest data by corresponding covariates
       group_by(across(all_of(c(strat_vars, covariate_var, "age")))) |>
       group_nest(keep = TRUE, .key = "data_fit") |>
       # create a nested list of fitting values
       mutate(data_fit = map(data_fit, ~ .x |>
                               dplyr::select(-c(all_of(strat_vars))) |>
                               as.list())
              )
   }
   
   # ------------------------------------------------------------------- #
   # PT2 - Estimating the MSM model and prepare the basis for prediction data
   # create fit model with arbitrary covariates
   # always assume that age splines are used and age is provided
   # NOTE: Takes time
   base_covs   <- names(spline_basis)
   all_covs    <- c(base_covs, covariate_var)
   cov_formula <- reformulate(all_covs)
   
   # calculate model
   result <- result |>
     group_by(across(all_of(strat_vars))) |>
     group_nest() |>
     # just in case
     na.omit() |>
     # removing solitary observations
     mutate(data = map(data, ~ .x |> 
                         group_by(hhidpn) |> 
                         filter(n() > 1) |> 
                         ungroup())) |>
     # modeling
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
   
   # prediction data basis
   pred <- prediction_grid |>
     # add our age interval
     mutate(age_interval = age_int) |>
     # join the model
     left_join(result, by = strat_vars)
   
   # ------------------------------------------------------------------- #
   # PT3 - This part is for CI calculation. assign variables to global env
   # If the method is normal or boot. msm assumes that all the variables
   # we have used in modelling and etc. are available in global .env.
   # this is due to windows parralelization routine
   # So I assign the variables to global and then remove them when we finish
   assign("cov_formula", cov_formula, envir = .GlobalEnv)
   assign("Q",           Q,           envir = .GlobalEnv)
   assign("ci",          ci_type,     envir = .GlobalEnv)
   assign("cl",          conf_level,  envir = .GlobalEnv)
   assign("B",           B,           envir = .GlobalEnv)
   assign("cores",       n_cores,     envir = .GlobalEnv)
 
   # ------------------------------------------------------------------- #
   # PT4 Calculate confidence intervals
   # There are 3 options now none, normal and bootstrapped
   # if none, we can calculate them directly and very fast using msm functions
   # Global assignment is not used here
   if(ci_type == "none") {
     
     result_df <- pred |>
       # First we calculate rate
       mutate(
         rate = map2(
           .x = model,
           .y = data_fit,
           ~ qmatrix.msm(x          = .x, 
                         covariates = .y)$estimates |>
             as.table() |>
             as.data.frame() |>
             rename(from     = Var1, 
                    to       = Var2, 
                    estimate = Freq) |>
             mutate(type = "q")
         ),
         # Here we calculate probability
         prob = pmap(
           list(model, data_fit, age_interval),
           ~ pmatrix.msm(x          = ..1, 
                         t          = ..3, 
                         covariates = ..2) |>
             as.table() |>
             as.data.frame() |>
             rename(from     = Var1, 
                    to       = Var2, 
                    estimate = Freq) |>
             mutate(type = "p")
         )) |>
       mutate(combined_tidy = map2(rate, prob, bind_rows)) |>
       select(-c(data_fit, age_interval, data, model, rate, prob)) |>
       unnest(combined_tidy) |>
       # remove D-D transitions (empty)
       filter(from != "State 3") |>
       # remove recovery possibility (empty)
       filter(!(from == "State 2" & to == "State 1"))
     
   } 
   
   if(ci_type == "normal") {
     
     # here the global assignment is in play. Except for cores argument
     # But the speed lost so so small, that it is best to follow the unified way
     pred$q_list <- map2(
       pred$model,
       pred$data_fit,
       ~ qmatrix.msm(x          = .x, 
                     covariates = .y,
                     ci         = ci,
                     cl         = cl,
                     # no cores argument here
                     B          = B))
     
     # this part is nothing but data cleaning and wrangling
     # used to return a dataframe with L estimate and U columns from q_list
     result_df <- pred |>
       mutate(
         # this part for q matrix
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
             left_join(upper,    by = c("from", "to")) |>
             mutate(type = "q")
         }
         ))  |> 
       mutate(
         # here we calculate p matrix from q using matrix exponentiation
         p_tidy = map2(q_list, age_interval, ~ {
           lower <- expm(.x$L * .y) |>
             as.data.frame() |>
             rownames_to_column("from") |>
             pivot_longer(-from, names_to = "to", values_to = "lower")
           
           estimate <- expm(.x$estimate * .y) |>
             as.data.frame() |>
             rownames_to_column("from") |>
             pivot_longer(-from, names_to = "to", values_to = "estimate")
           
           upper <- expm(.x$U * .y) |>
             as.data.frame() |>
             rownames_to_column("from") |>
             pivot_longer(-from, names_to = "to", values_to = "upper")
           
           lower |>
             left_join(estimate, by = c("from", "to")) |>
             left_join(upper, by = c("from", "to")) |>
             mutate(type = "p")
         }))  |>
       # bind both dataframes
       mutate(
         combined_tidy = map2(q_tidy, p_tidy, bind_rows)) |>
       # remove intermediate columns
       select(-c(data_fit, age_interval, data, 
                 model, q_list, q_tidy, p_tidy)) |> 
       # unnest
       unnest(combined_tidy) |>
       # remove D-D transitions (empty)
       filter(from != "State 3") |>
       # remove recovery possibility (empty)
       filter(!(from == "State 2" & to == "State 1")) 
     
   } 
   
   if(ci_type == "bootstrap") {
     # The most complex part that uses boots.
     # here cores arguments is used
     # takes an ton of time to calculate
     # here we calculate the rate list
     # This is the only part that uses bootstrapping and global variables
     # NOTE: Takes a very long time!
     
     # separate nested df into 2 lists
     # This supposedly facilitates the bootstrap fitting
     model    <- pred$model
     data_fit <- pred$data_fit
     
     # calculate CI
     q_list <- map2(
       model,
       data_fit,
       ~ qmatrix.msm(x          = .x,
                     covariates = .y,
                     ci         = ci,
                     cl         = cl,
                     B          = B,
                     cores      = cores))
     
     # paste q_list back as a corresponding column
     pred$q_list <- q_list
     
     # the rest is the same as in normal CI
     result_df <- pred |>
       mutate(
         # clean and prepare q data
         q_tidy = map(q_list, ~ {
           lower <- .x$L |>
             as.data.frame() |>
             rownames_to_column("from") |>
             pivot_longer(-from, names_to = "to", values_to = "lower")
           
           estimate <- .x$estimate |>
             as.data.frame() |>
             rownames_to_column("from") |>
             pivot_longer(-from, names_to = "to", values_to = "estimate")
           
           upper <- .x$U |>
             as.data.frame() |>
             rownames_to_column("from") |>
             pivot_longer(-from, names_to = "to", values_to = "upper")
           
           lower |>
             left_join(estimate, by = c("from", "to")) |>
             left_join(upper, by = c("from", "to")) |>
             mutate(type = "q")
         }),
         # use map2 for row-wise age_interval for p calculation
         p_tidy = map2(q_list, age_interval, ~ {
           lower <- expm(.x$L * .y) |>
             as.data.frame() |>
             rownames_to_column("from") |>
             pivot_longer(-from, names_to = "to", values_to = "lower")
           
           estimate <- expm(.x$estimate * .y) |>
             as.data.frame() |>
             rownames_to_column("from") |>
             pivot_longer(-from, names_to = "to", values_to = "estimate")
           
           upper <- expm(.x$U * .y) |>
             as.data.frame() |>
             rownames_to_column("from") |>
             pivot_longer(-from, names_to = "to", values_to = "upper")
           
           lower |>
             left_join(estimate, by = c("from", "to")) |>
             left_join(upper, by = c("from", "to")) |>
             mutate(type = "p")
         })) |>
       # bind
       mutate(
         combined_tidy = map2(q_tidy, p_tidy, bind_rows)
       ) |>
       # remove helper columns
       select(-c(data_fit, age_interval, data, model, 
                 q_list, q_tidy, p_tidy)) |>
       # unnest
       unnest(combined_tidy) |>
       # remove unnecessary transitions
       filter(from != "State 3") |>
       filter(!(from == "State 2" & to == "State 1"))
     
   }
   # ------------------------------------------------------------------- #
   # remove temporary variables from global env
   on.exit(rm("cov_formula",   envir = .GlobalEnv), add = TRUE)
   on.exit(rm("Q",             envir = .GlobalEnv), add = TRUE)
   on.exit(rm("ci",            envir = .GlobalEnv), add = TRUE)
   on.exit(rm("cl",            envir = .GlobalEnv), add = TRUE)
   on.exit(rm("B",             envir = .GlobalEnv), add = TRUE)
   on.exit(rm("cores",         envir = .GlobalEnv), add = TRUE)
   
   # return the final dataframe
   return(result_df)
   
 }
# ------------------------------------------------------------------- #
# second version working in new wrapper...
# breaks!
fit_msm2 <- function(.data,
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
  # --- PT3 + PT4 MODIFIED -------------------------------------------- #
  if (ci_type == "none") {

    pred <- prediction_grid |>
      mutate(age_interval = age_int) |>
      left_join(result, by = strat_vars) |>
      mutate(qp = map2(
        .x = model,
        .y = data_fit,
        ~ qmatrix.msm_wrapper(x = .x,
                              covariates = .y,
                              ci = "none",
                              age_interval = age_int)
      )) |>
      mutate(
        rate = map(qp, ~ .x$Q$estimate |>
                     as.table() |>
                     as.data.frame() |>
                     rename(from = Var1, to = Var2, rate = Freq)),
        prob = map(qp, ~ .x$P$estimate |>
                     as.table() |>
                     as.data.frame() |>
                     rename(from = Var1, to = Var2, prob = Freq))
      )

    rate <- pred |>
      dplyr::select(all_of(c(strat_vars, covariate_var)), age, rate) |>
      unnest(rate)

    prob <- pred |>
      dplyr::select(all_of(c(strat_vars, covariate_var)), age, prob) |>
      unnest(prob)

    result_df <- rate |>
      full_join(prob) |>
      filter(from != "State 3") |>
      filter(!(from == "State 2" & to == "State 1"))

  } else {
    pred <- prediction_grid |>
      mutate(age_interval = age_int) |>
      left_join(result, by = strat_vars) |>
      mutate(ci    = ci_type,
             cl    = conf_level,
             B     = B,
             cores = n_cores) |>
      mutate(Q = list(Q)) |>  # pass Q matrix just in case needed inside wrapper
      mutate(qp_list = pmap(
        list(model, data_fit, ci, cl, B, cores),
        ~ qmatrix.msm_wrapper(
          x            = ..1,
          covariates   = ..2,
          ci           = ..3,
          B            = ..5,
          cores        = ..6,
          age_interval = age_int
        )
      ))

    result_df <- pred |>
      mutate(
        q_tidy = map(qp_list, ~ {
          lower <- .x$Q$L |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from, names_to = "to", values_to = "lower")

          estimate <- .x$Q$estimate |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from, names_to = "to", values_to = "estimate")

          upper <- .x$Q$U |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from, names_to = "to", values_to = "upper")

          lower |>
            left_join(estimate, by = c("from", "to")) |>
            left_join(upper, by = c("from", "to")) |>
            mutate(type = "q")
        }),
        p_tidy = map(qp_list, ~ {
          lower <- .x$P$L |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from, names_to = "to", values_to = "lower") |>
            mutate(from = str_c("State", from, sep = " "),
                   to   = str_c("State", parse_number(to), sep = " "))

          estimate <- .x$P$estimate |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from, names_to = "to", values_to = "estimate")

          upper <- .x$P$U |>
            as.data.frame() |>
            rownames_to_column("from") |>
            pivot_longer(-from, names_to = "to", values_to = "upper") |>
            mutate(from = str_c("State", from, sep = " "),
                   to   = str_c("State", parse_number(to), sep = " "))

          lower |>
            left_join(estimate, by = c("from", "to")) |>
            left_join(upper, by = c("from", "to")) |>
            mutate(type = "p")
        }),
        combined_tidy = map2(q_tidy, p_tidy, bind_rows)
      ) |>
      select(-c(data_fit, age_interval, data, model,
                qp_list, q_tidy, p_tidy, ci, cl, B, cores, Q)) |>
      unnest(combined_tidy) |>
      filter(from != "State 3") |>
      filter(!(from == "State 2" & to == "State 1"))
  }
  
  # return final result
  return(result_df)
}

# --------------------------------------------------------------- # 
# new function for parallel calculation

fit_msm <- function(.data, # .data for piping
                    # main variables
                    # independent stratas
                    strat_vars    = NULL,
                    # covariates in the msm model
                    covariate_var = NULL,
                    # only used for continuous variables
                    cont_grid     = NULL,
                    # age grid specification part from, to
                    age_from_to   = c(50, 100),
                    # age increment
                    age_int       = 0.25,
                    # spline specification part
                    # do we want to calculate spline at all?
                    calc_spline   = TRUE,
                    # what type of spline? ns or bs
                    spline_type   = "ns",
                    # number of splines
                    spline_df     = NULL,
                    # conf.int calculation part
                    # do we want to calculate CI at all?
                    ci            = FALSE,
                    # if yes, then what number of parralel processes
                    # do we want to use for msm modelling?
                    n_cores       = 10,
                    # if yes specify the number of bootstraps we
                    # want to calculate our CI with
                    B             = 2,
                    # specify confidence level here
                    conf_level    = 0.95,
                    # initialize Q matrix
                    Q = rbind(
                      c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
                      c(0, -.01,   0.1),  # dementia can go to death
                      c(0, 0,   0)        # death is absorbing
                    )) {
  
  # ------------------------------------------------------------------- #
  # PT1 - Creating the FIT and PREDICT DATASETS
  # Allows splines and no splines.
  # Allows the use of continuous covariate like year e.g. (2010, 2012, 2020)
  
  if(!is.null(cont_grid)) {
    
    vars <- map(c(strat_vars), ~ unique(na.omit(.data[[.x]])))
    vars[[(length(vars) + 1)]] <- cont_grid
    names(vars) <- c(strat_vars, covariate_var)
    
    # otherwise we use this part
  } else { 
    
    vars <- map(c(covariate_var, strat_vars),~ unique(na.omit(.data[[.x]]))) |>
      set_names(c(covariate_var, strat_vars))
    
  }
  
  # add base predictors. in our model it is age
  # from age_pred_grid from fisr to last by age_int
  base_grid <- list(
    age = seq(age_from_to[1], age_from_to[2], by = age_int)
  )
  
  # construct prediction grid dynamically
  prediction_grid <- exec(crossing, !!!c(base_grid, vars))
  
  if(calc_spline) { 
    
    # construct spline basis for age
    # Here we can also add another types of splines in future
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
    
    # bind data together
    prediction_grid <- prediction_grid |>
      # nest data by corresponding covariates
      group_by(across(all_of(c(strat_vars, covariate_var, "age")))) |>
      group_nest(keep = TRUE, .key = "data_fit") |>
      # create a nested list of fitting values
      mutate(data_fit = map(data_fit, ~ .x |>
                              dplyr::select(-c(all_of(strat_vars), age)) |>
                              as.list())
      )
    
  } else { # if we do not want to use splines
    # The difference here is in the data_fit. Previous removes age, this keeps it
    
    spline_basis <- prediction_grid |>
      dplyr::select(age)
    
    # here we assign original data for compatibility with other method
    result <- .data
    
    # bind data together
    prediction_grid <- prediction_grid |>
      # nest data by corresponding covariates
      group_by(across(all_of(c(strat_vars, covariate_var, "age")))) |>
      group_nest(keep = TRUE, .key = "data_fit") |>
      # create a nested list of fitting values
      mutate(data_fit = map(data_fit, ~ .x |>
                              dplyr::select(-c(all_of(strat_vars))) |>
                              as.list())
      )
  }
  
  # ------------------------------------------------------------------- #
  # PT2 - Specify weather we want ot calculate CI
  if(ci) { 
    
    # if we want, then we create bootstrap resamples with the prespecified B
    hrs_boot <- group_bootstraps2(result,
                                  group  = "hhidpn",
                                  weight = "wtcrnh",
                                  times  = B)
    
    # we collect these resamples into single data
    hrs_boot <- map(hrs_boot$splits, ~ .x |>
                      analysis() |>
                      # Ensure ordering by subject and time
                      arrange(hhidpn, age)
    )
    
    # and bind them together creating boot id (for parralel)
    hrs_boot <- hrs_boot |>
      bind_rows(.id = "boot")
    
  } else { 
    
    # in other case we have just one boot and we assign B to 1
    # I do this for code unification
    hrs_boot <- result |>
      mutate(boot = 1)
    
    B <- 1
    
  }
  
  # ------------------------------------------------------------------- #
  # PT3 - Estimating the MSM model and prepare the basis for prediction data
  # create fit model with arbitrary covariates
  # always assume that age splines are used and age is provided
  # NOTE: Takes time
  
  # this first part if being calculated anyway
  # it creates the predict data essentially
  base_covs   <- names(spline_basis)
  all_covs    <- c(base_covs, covariate_var)
  cov_formula <- reformulate(all_covs)
  
  # calculate model
  result <- hrs_boot |>
    group_by(across(all_of(c("boot", strat_vars)))) |>
    group_nest() |>
    # just in case
    na.omit() |>
    # removing solitary observations
    mutate(data = map(data, ~ .x |> 
                        group_by(hhidpn) |> 
                        filter(n() > 1) |> 
                        ungroup())) |>
    # these steps are used for future parallel processing
    unnest(data)
  
  
  # here we use parrallel processing to fir our msm models
  # I use multidplyr currently
  if(ci) {
    
    # make cluster
    cluster <- new_cluster(n_cores)
    # export packages to cluster
    cluster_library(cluster, c("msm", "dplyr", "tidyselect"))
    
    # export necessary definitions to each cluster
    cluster_assign(cluster, Q = Q)
    cluster_assign(cluster, cov_formula = cov_formula)
    cluster_assign(cluster, strat_vars  = strat_vars)
    
    # partition data into clusters by intersection of boot number and strat_vars
    result <- result |>
      group_by(across(all_of(c("boot", strat_vars)))) |>
      partition(cluster)
    
    # fit the model in parallel
    models <- result |>
      do(mod = msm(
        formula    = state_msm ~ age,
        subject    = hhidpn,
        data       = .data,
        qmatrix    = Q,
        obstype    = obstype,
        deathexact = 3,
        covariates = cov_formula,
        control    = list(fnscale = 5000, 
                          maxit   = 25000), # ?????
        gen.inits  = TRUE,
        method     = "BFGS"
      )) |> 
      # this one is for binding
      mutate(boot = as.numeric(boot))|>
      # collect data
      collect()
    
  } else { 
    
    models <- result |>
      group_by(across(all_of(c("boot", strat_vars)))) |>
      do(mod = msm(
        formula    = state_msm ~ age,
        subject    = hhidpn,
        data       = .data,
        qmatrix    = Q,
        obstype    = obstype,
        deathexact = 3,
        covariates = cov_formula,
        control    = list(fnscale = 5000, 
                          maxit   = 25000), #?????
        gen.inits  = TRUE,
        method     = "BFGS"
      )) |>
      mutate(boot = as.numeric(boot))
    
  }
  
  # prediction data basis
  # create copies of data for each bootstrap for CI
  pred <- crossing(prediction_grid, boot = c(1:B)) |>
    # join with the model
    left_join(models, by = c("boot", strat_vars)) |>
    # calculate q_matrices
    group_by(across(all_of(c(strat_vars)))) |>
    mutate(q_list = map2(mod, data_fit, ~
                           qmatrix.msm(x          = .x,
                                       covariates = .y)$estimate))
  
  # ------------------------------------------------------------------- #
  # PT4 - Caculate CI.
  # if true
  if(ci) {
    
    # define upper and lower prob
    alpha      <- 1 - conf_level
    lower_prob <- alpha     / 2      
    upper_prob <- 1 - alpha / 2
    
    # tidy q and calculate p
    result_df <- pred |>
      # for each age covariate and strata
      group_by(across(all_of(c(strat_vars, covariate_var, "age")))) |>
      # calculate q array for all bootstraps
      summarise(
        # Collect all bootstrap matrices into a 3D array
        q_array = list(simplify2array(q_list)), .groups = "drop"
      ) |>
      # here we calculate the estimate L and U from for Q our q rray
      mutate(
        q_mean  = map(q_array, ~ apply(.x, c(1, 2), mean)),
        # note use of upper and lower
        q_lower = map(q_array, ~ apply(.x, c(1, 2), quantile, probs = lower_prob)),
        q_upper = map(q_array, ~ apply(.x, c(1, 2), quantile, probs = upper_prob))
      ) |>
      # this part here is just making q matrices a tidy df
      # basically the same for estimate L and U
      mutate(
        mean_q = map(q_mean, ~ .x |>
                       as.table() |> 
                       as.data.frame() |> 
                       rename(from     = Var1,
                              to       = Var2,
                              estimate = Freq)),
        lower_q = map(q_lower, ~ .x |> 
                        as.table() |> 
                        as.data.frame() |> 
                        rename(from  = Var1, 
                               to    = Var2, 
                               lower = Freq)),
        upper_q = map(q_upper, ~ .x |>
                        as.table() |> 
                        as.data.frame() |> 
                        rename(from  = Var1, 
                               to    = Var2, 
                               upper = Freq)),
        # Combine all q data into a single data frame
        combined_q = pmap(
          list(mean_q, lower_q, upper_q), 
          ~ reduce(list(..1, ..2, ..3), left_join, by = c("from", "to")))
      ) |>
      # now we calclate p matrices with expm and our age_int
      mutate(
        # exponentiate q matrix and assign age_int
        mean_p  = map(q_mean, ~ expm(.x * age_int) |>
                        as.table() |>
                        as.data.frame() |>
                        rename(from       = Var1,
                               to         = Var2,
                               p_estimate = Freq)),
        # same for L and U
        lower_p = map(q_lower, ~ expm(.x * age_int) |>
                        as.table() |>
                        as.data.frame() |>
                        rename(from    = Var1, 
                               to      = Var2, 
                               p_lower = Freq)),
        upper_p = map(q_upper, ~ expm(.x * age_int) |>
                        as.table() |>
                        as.data.frame() |>
                        rename(from    = Var1, 
                               to      = Var2, 
                               p_upper = Freq)),
        # combine all p into single df
        combined_p = pmap(
          list(mean_p, lower_p, upper_p), 
          ~ reduce(list(..1, ..2, ..3), left_join, by = c("from", "to")))
      ) |>
      # combine p and q into one df
      mutate(combined = map2(combined_q, combined_p, ~ .x |>
                               full_join(.y, by = c("from", "to")))
      ) |>
      # select only columns that we need
      select(all_of(c(strat_vars, covariate_var, "age", "combined"))) |>
      # unnest combined df
      unnest(combined)|>
      # remove D-D transitions (empty)
      filter(from != "State 3") |>
      # remove recovery possibility (empty)
      filter(!(from == "State 2" & to == "State 1"))
    
    # if no CI needed
  } else {
    
    # no CI process is simpler
    result_df <- pred |>
      # First we calculate q list
      mutate(
        q_list = map2(
          .x = mod,
          .y = data_fit,
          ~ qmatrix.msm(x          = .x, 
                        covariates = .y)$estimates
        )) |>
      # tidy it and make it a tibble
      mutate(mean_q = map(q_list, ~ .x |>
                            as.table() |>
                            as.data.frame() |>
                            rename(from       = Var1, 
                                   to         = Var2,
                                   q_estimate = Freq))) |>
      # calculate p with expm from q and tidy it
      mutate(
        mean_p  = map(q_list, ~ expm(.x * age_int) |>
                        as.table() |>
                        as.data.frame() |>
                        rename(
                          from       = Var1,
                          to         = Var2,
                          p_estimate = Freq
                        ))) |>
      # here we combine p and q tibbles into one df
      mutate(combined = map2(mean_q, mean_p, ~ .x |>
                               left_join(.y, by = c("from", "to")))) |>
      # remove helper columns
      select(-c(data_fit, mod, mean_q, mean_p, q_list, boot)) |>
      # unnest combined df
      unnest(combined) |>
      # remove D-D transitions (empty)
      filter(from != "State 3") |>
      # remove recovery possibility (empty)
      filter(!(from == "State 2" & to == "State 1"))
    
  }
  
  # return the final dataframe
  return(result_df)
  
}
