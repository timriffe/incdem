library(rsample)
library(tidvyerse)
# maybe pre-define which ids enter each boot replicate, using a minimal
# data structure
source("R/01_prepare_hrs.R")
source("R/zzz_boot_test_functions.R")

hrs_boot <- group_bootstraps2(hrs_to_fit,
                              group = "hhidpn",
                              weight = "wtcrnh",
                              times = 10)

# prepare parralelization
set.seed(123)
n_cores <- 7
registerDoParallel(cores = n_cores)

#  Q matrix
Q <- rbind(
  c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
  c(0,   -.01, 0.1),  # dementia can go to death
  c(0,    0,   0)     # death is absorbing
)

# Prepare  list of bootstrap samples
example <- map(hrs_boot$splits, analysis)


# Function to fit MSM model on a single dataset
fit_on_split <- function(split_df) {
  # here is it
  split_df <- split_df %>%
    arrange(hhidpn, age)  # Ensure ordering by subject and time
  
  msm_model <- msm(
    state_msm ~ age,
    subject   = hhidpn,
    data      = split_df,
    qmatrix   = Q,
    covariates = ~ female,
    control   = list(fnscale = 5000)
  )
  
  qmatrix.msm(msm_model)
}

# Run in parallel
results <- foreach(split_df = example,
                   .packages = c("msm", "dplyr")) %dopar% {
                     fit_on_split(split_df)
                   }

stopImplicitCluster()

# collect
q_list  <- map(results, ~ .x$estimate)
q_array <- simplify2array(q_list)

# Calculate element-wise quantiles
q_lower <- apply(q_array, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
q_upper <- apply(q_array, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
q_mean  <- apply(q_array, c(1, 2), mean, na.rm = TRUE)  # Optional

# Combine into a named list
q_ci <- list(
  estimate = q_mean,
  L = q_lower,
  U = q_upper
)
