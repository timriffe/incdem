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

# ------------------------------------------------------------------- #
source("R/01_data_preparation.R")
source("R/02_estimate_models.R")
# ------------------------------------------------------------------- #
# set desired age interval
age_int <- 0.25
# construct prediction grid
prediction_grid <- crossing(
  sex = c("Male", "Female"),
  age = seq(50, 100, by = age_int),
  obs_date = c(2000, 2010, 2020)
)
spline_basis <- predict(spline_basis_fit, 
                        newx = prediction_grid$age) |>
  as.data.frame() |>
  setNames(paste0("age_spline", 1:3))
prediction_grid <- bind_cols(prediction_grid, spline_basis)
# ------------------------------------------------------------------- #
# rename the model sex for binding
model <- split_data |>
  mutate(female = ifelse(female == 1, "Female", "Male"))

# this is the full fit for prob and rates
pred <- prediction_grid |>
  # for each combination of these variables
  group_nest(sex, age, obs_date, .key = "data_fit") |>
  # create the fitting list
  mutate(data_fit = map2(
    .x = data_fit,
    .y = obs_date,
    ~ mutate(.x, obs_date = .y) |>
      as.list()
  )) |>
  # join the model
  left_join(model, by = c("sex" = "female")) |>
  # set age intervals
  mutate(age_interval = age_int) |> # set age interval here
  # create rate
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
    # create probability
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
# unnest rate
rate <- pred |>
  dplyr::select(sex, age, obs_date, rate) |>
  unnest(rate)

# unnest probability
prob <- pred |>
  dplyr::select(sex, age, obs_date, prob) |>
  unnest(prob)

# this is the final result
result_df <- rate |>
  full_join(prob, 
            by = join_by(sex, age,
                         obs_date, from, to)) |>
  # remove D-D transitions (empty)
  filter(from != "State 3") |>
  # remove recovery possibility (empty)
  filter(!(from == "State 2" & to == "State 1"))
# ------------------------------------------------------------------- #
# write to csv.gz and not Rdata
write_csv(result_df, file = "Data/results.csv.gz")
