# ------------------------------------------------------------------- #
source("R/01b_data_preparation.R")
# ------------------------------------------------------------------- #
# Q tells msm what the valid transitions are
Q <- rbind(
  c(0, 0.1, 0.1),  # healthy can go to dementia or death
  c(0, 0,   0.1),  # dementia can go to death
  c(0, 0,   0)     # death is absorbing
)
diag(Q) <- -rowSums(Q)
# ------------------------------------------------------------------- #
# fit separate models to males and females
# since there is no by argument in msm (TRUE)
# lets use newly created obs_date. We will have more observations
split_data <- hrs_msm |>
  # a test condition to check pandemic effect on mort trend
  mutate(birth_date         = as_date(birth_date, origin = "1960-01-01"),
         birth_date_decimal = decimal_date(birth_date),
         obs_date           = birth_date_decimal + age) |> 
  filter(obs_date < (as_date("2019-jan-01") |> decimal_date())) |>
  #filter(obs_date < (as_date("2019-dec-31") |> decimal_date())) |> 
  # filter(obs_date < (as_date("2020-feb-28") |> decimal_date())) |> 
  group_by(female) |>
  group_nest() |>
  mutate(model = map(data, ~ msm(
    state_msm ~ age,
    subject    = hhidpn,
    data       = .x,
    qmatrix    = Q,
    obstype    = obstype,
    deathexact = 3,
    covariates = ~ age_spline1 + age_spline2 + age_spline3 + obs_date,
    control    = list(fnscale = 5000, maxit = 25000),
    gen.inits  = TRUE,
    method     = "BFGS" 
  )))
# ------------------------------------------------------------------- #
# Checking the model
# 0 means it converged
split_data$model[[1]]$opt$convergence
split_data$model[[2]]$opt$convergence

# check transition intensity matrices
# H-D and U-D (should be larger) and both > 0
qmatrix.msm(split_data$model[[1]])
qmatrix.msm(split_data$model[[2]])

# probability matrix
pmatrix.msm(split_data$model[[1]], t = 1)
pmatrix.msm(split_data$model[[2]], t = 1)

# Looks like female e(50) too high?
# But not sure how it's calculated;
# we should compare with our own calculations to be sure.
split_data$model[[1]]$sojourn
split_data$model[[2]]$sojourn

# extract observed and fitted prevalence matrices
# we can plot them later if needed
# NOTE: Takes long time
# prevalence.msm(split_data$model[[1]])
# prevalence.msm(split_data$model[[2]])

# expected time in each state
efpt.msm(split_data$model[[1]], tostate = 3, covariates = "mean")
efpt.msm(split_data$model[[2]], tostate = 3, covariates = "mean")
# ------------------------------------------------------------------- #