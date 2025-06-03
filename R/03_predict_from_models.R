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
  int_date_decimal = c(2000, 2010, 2020)
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
  group_nest(sex, age, int_date_decimal, .key = "data_fit") |>
  # create the fitting list
  mutate(data_fit = map2(
    .x = data_fit,
    .y = int_date_decimal,
    ~ mutate(.x, int_date_decimal = .y) |>
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
  dplyr::select(sex, age, int_date_decimal, rate) |>
  unnest(rate)

# unnest probability
prob <- pred |>
  dplyr::select(sex, age, int_date_decimal, prob) |>
  unnest(prob)

# this is the final result
result_df <- rate |>
    full_join(prob, 
              by = join_by(sex, age,
                           int_date_decimal, from, to)) |>
  # remove D-D transitions (empty)
  filter(from != "State 3") |>
  # remove recovery possibility (empty)
  filter(!(from == "State 2" & to == "State 1"))
# ------------------------------------------------------------------- #
# write to csv.gz and not Rdata
write_csv(result_df, file = "Data/results.csv.gz")
