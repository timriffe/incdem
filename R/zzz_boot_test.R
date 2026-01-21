# library(rsample)
# library(tidvyerse)
# # maybe pre-define which ids enter each boot replicate, using a minimal
# # data structure
# source("R/01_prepare_hrs.R")
# source("R/zzz_boot_test_functions.R")
# 
# hrs_boot <- group_bootstraps2(hrs_to_fit,
#                               group = "hhidpn",
#                               weight = "wtcrnh",
#                               times = 10)
# 
# # prepare parralelization
# set.seed(123)
# n_cores <- 7
# registerDoParallel(cores = n_cores)
# 
# #  Q matrix
# Q <- rbind(
#   c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
#   c(0,   -.01, 0.1),  # dementia can go to death
#   c(0,    0,   0)     # death is absorbing
# )
# 
# # Prepare  list of bootstrap samples
# example <- map(hrs_boot$splits, analysis)
# 
# 
# # Function to fit MSM model on a single dataset
# fit_on_split <- function(split_df) {
#   # here is it
#   split_df <- split_df %>%
#     arrange(hhidpn, age)  # Ensure ordering by subject and time
#   
#   msm_model <- msm(
#     state_msm ~ age,
#     subject   = hhidpn,
#     data      = split_df,
#     qmatrix   = Q,
#     covariates = ~ female,
#     control   = list(fnscale = 5000)
#   )
#   
#   qmatrix.msm(msm_model)
# }
# 
# # Run in parallel
# results <- foreach(split_df = example,
#                    .packages = c("msm", "dplyr")) %dopar% {
#                      fit_on_split(split_df)
#                    }
# 
# stopImplicitCluster()
# 
# # collect
# q_list  <- map(results, ~ .x$estimate)
# q_array <- simplify2array(q_list)
# 
# # Calculate element-wise quantiles
# q_lower <- apply(q_array, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
# q_upper <- apply(q_array, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
# q_mean  <- apply(q_array, c(1, 2), mean, na.rm = TRUE)  # Optional
# 
# # Combine into a named list
# q_ci <- list(
#   estimate = q_mean,
#   L = q_lower,
#   U = q_upper
# )
# --------------------------------------------------------------------------- # 
# here I show the parallel msm results
# functions
source("R/01_prepare_hrs.R")
source("R/zzz_boot_test_functions.R")

# data with columns that we use
hrs_to_fit_short <- hrs_to_fit |>
  mutate(period = case_when(between(obs_date, 2000, 2006) ~ "period 1",
                            between(obs_date, 2006, 2012) ~ "period 2",
                            between(obs_date, 2012, 2020) ~ "period 3")) |>
  select(hhidpn, female, education,
         wtcrnh, age, hypertension:stroke,
         state_msm, ever_dementia, obstype, period)

# NO CI
no_boot <- hrs_to_fit_short |>
  fit_msm(strat_vars    = c("female"), 
          covariate_var = c("period"), 
          age_from_to   = c(50, 100), 
          age_int       = 0.25,
          cont_grid     = NULL,
          calc_spline   = TRUE, # if false it does not matter what spline_df and spline_type are
          spline_df     = 3,    # only for calc_spline = T
          spline_type   = "ns", # only for calc_spline = T
          # CI part
          # no ci in this case
          ci            = FALSE, # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          # ignored if  ci = F
          n_cores       = 8, 
          # ignored if ci = F
          B             = 10,
          # ignored  if ci = F
          conf_level    = 0.95,
          Q = rbind(
            c(-0.2, 0.1, 0.1),
            c(0, -.01,   0.1), 
            c(0, 0,   0) 
          ))

# WITH CI
with_boot <- hrs_to_fit_short |>
  fit_msm(strat_vars    = c("female"), 
          covariate_var = c("period"),
          age_from_to   = c(50, 100), 
          age_int       = 0.25,
          cont_grid     = NULL, 
          calc_spline   = TRUE,
          spline_df     = 3, 
          spline_type   = "ns", 
          # here we calculate CI
          ci            = TRUE, # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          # next 3 arguments are used
          n_cores       = 5, 
          # notice it is very low for testing
          B             = 30, 
          conf_level    = 0.95,
          Q = rbind(
            c(-0.2, 0.1, 0.1), 
            c(0, -.01,   0.1),  
            c(0, 0,   0)
          ))
