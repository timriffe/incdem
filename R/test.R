# install.packages("msm")
library(msm)
library(tidyverse)
library(lubridate)
library(splines)
hrs <- read_csv("Data/riffe_incdem_20250522.csv")



# first pass processing
hrs_msm <-
  hrs |> 
  # pick age range to fit to, based on plot of support.
  # if all ages 
  filter(between(age,55,97)) |> 
  mutate(
    int_date = suppressWarnings(as.integer(int_date)),
    interview_date = lubridate::as_date(int_date, origin = "1960-01-01"),
    int_date_decimal = decimal_date(interview_date)
  ) |> 
  arrange(hhidpn, age) |> 
  # mutate(age_diff = age - (int_date_decimal - decimal_date(as_date(birth_date, origin = "1960-01-01")))) |> 
  mutate(
    state_msm = state + 1  # msm expects states starting at 1
  ) |> 
  group_by(hhidpn) |> 
  filter(n() > 1) |> # supposed to remove solitary observations.
  ungroup()

# all filtering before here; needed to create spline 
# basis with consistent age range in fit and in predictions.
# spline_basis_fit is also used later for this reason
spline_basis_fit   <- ns(hrs_msm$age, df = 3)
age_splines        <- as.data.frame(spline_basis_fit)
names(age_splines) <- paste0("age_spline", 1:3)

# This second part is supposed to eliminate impossible transitions.
hrs_msm <-
  hrs_msm |> 
  bind_cols(age_splines) |> 
  arrange(hhidpn, age) |> group_by(hhidpn) |> 
  arrange(age) |> 
  mutate(
    ever_dementia = cumany(state == 1),
    state_clean = case_when(
      state == 2 ~ 2,  # preserve death
      ever_dementia & state == 0 ~ 1,  # impute "recovery" as still dementia
      TRUE ~ state
    )
  ) |> 
  ungroup() |> 
  group_by(hhidpn)  |> 
  arrange(age) |> 
  mutate(
    died = cumany(state_clean == 2)
  ) |>
  filter(!(died & state_clean != 2)) |> 
  ungroup() |> 
  arrange(hhidpn,age) |> 
  # treat death times as exact, but othe transition
  # times as unknown.
  mutate(obstype = ifelse(state_msm == 3,3,1))

hrs_msm <- hrs_msm |>
  group_by(hhidpn) |> 
  # supposed to remove solitary observations.
  filter(n() > 1) |>
  ungroup()

# Q tells msm what the valid transitions are
Q <- rbind(
  c(0, 0.1, 0.1),  # healthy can go to dementia or death
  c(0, 0,   0.1),  # dementia can go to death
  c(0, 0,   0)     # death is absorbing
)
diag(Q) <- -rowSums(Q)

# fit separate models to males and females
# since there is no by argument in msm...
model_msm_female <- msm(
  state_msm ~ age,
  subject = hhidpn,
  data = filter(hrs_msm, female == 1),
  qmatrix = Q,
  obstype = obstype,
  deathexact = 3,
  covariates = ~ age_spline1 + age_spline2 + age_spline3 + int_date_decimal,
  control = list(fnscale = 5000, maxit = 20000),
  gen.inits = TRUE,
  method = "BFGS"
)

model_msm_male <- msm(
  state_msm ~ age,
  subject = hhidpn,
  data = filter(hrs_msm, female == 0),
  qmatrix = Q,
  obstype = obstype,
  deathexact = 3,
  covariates = ~ age_spline1 + age_spline2 + age_spline3 + int_date_decimal,
  control = list(fnscale = 5000, maxit = 20000),
  gen.inits = TRUE,
  method = "BFGS"
)

# 0 means it converged
model_msm_female$opt$convergence
model_msm_male$opt$convergence

# Looks like female e(50) too high?
# But not sure how it's calculated;
# we should compare with our own calcs to be sure.
model_msm_female$sojourn
model_msm_male$sojourn


# stick together to pass in
models <- list(Male = model_msm_male,
               Female = model_msm_female)

# For the prediction grid, we need to be careful to give
# splines matched to age in the same way as in the hrs
# data above. If we add covariates, then these need to be 
# accounted for here somehow.
prediction_grid <- crossing(
  sex = c("Male", "Female"),
  age = seq(50, 100, by = 0.25),
  int_date_decimal = c(2000, 2010, 2020)
)
spline_basis <- predict(spline_basis_fit, 
                        newx = prediction_grid$age) |>
  as.data.frame() |>
  setNames(paste0("age_spline", 1:3))
prediction_grid <- bind_cols(prediction_grid, spline_basis)

  
# custom function to predict haz and prob for our grid.
# this also lacks time still. Note we pass in models
# rather than counting on up-scoping.
#' @param .x a row of data excluding key variables
#' @param .y a row of key variables
#' @param age_interval same as time interval, default .25 (quarters)
#' @param models list with elements "Male" and "Female", the fitted msm models.
predict_hazard_and_prob <- function(.x, 
                                    .y, 
                                    age_interval = 0.25, 
                                    models) {
  model <- models[[.y$sex]]
  
  covariates_row <- .x |>
    select(starts_with("age_spline")) |>
    slice(1) |>
    as.list()
  names(covariates_row) <- c(paste0("age_spline", 1:3))
  covariates_row$int_date_decimal <- .y$int_date_decimal
  
  
  qmat <- qmatrix.msm(model, 
                      covariates = covariates_row)$estimates |> 
    as.table() |> 
    as.data.frame() |> 
    rename(from = Var1, to = Var2, rate = Freq)
  
  pmat <- pmatrix.msm(model, 
                      t = age_interval, 
                      covariates = covariates_row) |> 
    as.table() |> 
    as.data.frame() |> 
    rename(from = Var1, to = Var2, prob = Freq)
  
  left_join(
    qmat,
    pmat,
    by = join_by(from,to))
}

predict_hazard_and_prob <- function(.x, .y, age_interval = 0.25, models) {
  model <- models[[.y$sex]]
  
  covariates_row <- .x |>
    select(starts_with("age_spline")) |>
    slice(1) |>
    as.list()
  names(covariates_row) <- c(paste0("age_spline", 1:3))
  covariates_row$int_date_decimal <- .y$int_date_decimal
  
  # Hazards with CIs
  qmat_full <- qmatrix.msm(model, covariates = covariates_row, ci = "delta")
  qmat_est  <- as.data.frame(as.table(qmat_full$estimates))
  qmat_low  <- as.data.frame(as.table(qmat_full$L))
  qmat_high <- as.data.frame(as.table(qmat_full$U))
  
  q_df <- qmat_est |>
    rename(from = Var1, to = Var2, rate = Freq) |>
    mutate(lower_rate = qmat_low$Freq,
           upper_rate = qmat_high$Freq)
  
  # Probabilities with CIs
  pmat_full <- pmatrix.msm(model, t = age_interval,
                           covariates = covariates_row, ci = "normal")
  pmat_est  <- as.data.frame(as.table(pmat_full$estimates))
  pmat_low  <- as.data.frame(as.table(pmat_full$L))
  pmat_high <- as.data.frame(as.table(pmat_full$U))
  
  p_df <- pmat_est |>
    rename(from = Var1, to = Var2, prob = Freq) |>
    mutate(lower_prob = pmat_low$Freq,
           upper_prob = pmat_high$Freq)
  
  # Join both
  left_join(q_df, p_df, by = c("from", "to"))
}

# this one does normal cis for both rates and probs
predict_hazard_and_prob <- function(.x, .y, age_interval = 0.25, models) {
  model <- models[[.y$sex]]
  
  covariates_row <- .x |>
    select(starts_with("age_spline")) |>
    slice(1) |>
    as.list()
  names(covariates_row) <- c(paste0("age_spline", 1:3))
  covariates_row$int_date_decimal <- .y$int_date_decimal
  
  # Hazards with CIs
  qmat_full <- qmatrix.msm(model, covariates = covariates_row, ci = "normal")
  qmat_est  <- as.data.frame(as.table(qmat_full$estimates))
  qmat_low  <- as.data.frame(as.table(qmat_full$L))
  qmat_high <- as.data.frame(as.table(qmat_full$U))
  
  q_df <- qmat_est |>
    rename(from = Var1, to = Var2, rate = Freq) |>
    mutate(lower_rate = qmat_low$Freq,
           upper_rate = qmat_high$Freq)
  
  # Probabilities with CIs
  pmat_full <- pmatrix.msm(model, t = age_interval,
                           covariates = covariates_row, ci = "normal")
  pmat_est  <- as.data.frame(as.table(pmat_full$estimates))
  pmat_low  <- as.data.frame(as.table(pmat_full$L))
  pmat_high <- as.data.frame(as.table(pmat_full$U))
  
  p_df <- pmat_est |>
    rename(from = Var1, to = Var2, prob = Freq) |>
    mutate(lower_prob = pmat_low$Freq,
           upper_prob = pmat_high$Freq)
  
  # Join both
  left_join(q_df, p_df, by = c("from", "to"))
}

# this one does normal for rates, then uses those 3 matrices to get probs.
# second step not working...
predict_hazard_and_prob <- function(.x, .y, age_interval = 0.25, models) {
  model <- models[[.y$sex]]
  
  covariates_row <- .x |>
    select(starts_with("age_spline")) |>
    slice(1) |>
    as.list()
  names(covariates_row) <- c(paste0("age_spline", 1:3))
  covariates_row$int_date_decimal <- .y$int_date_decimal
  
  # Get estimated Q matrix and its CIs
  qmat_full <- qmatrix.msm(model, covariates = covariates_row, ci = "normal")
  qmat_est  <- qmat_full$estimates
  qmat_low  <- qmat_full$L
  qmat_high <- qmat_full$U
  
  # Predict transition probabilities at fixed Q matrices
  p_est <- pmatrix.msm(model, t = age_interval, qmatrix = qmat_est)
  p_low <- pmatrix.msm(model, t = age_interval, qmatrix = qmat_low)
  p_high <- pmatrix.msm(model, t = age_interval, qmatrix = qmat_high)
  
  # Format transition rates
  q_df <- as.data.frame(as.table(qmat_est)) |>
    rename(from = Var1, to = Var2, rate = Freq) |>
    mutate(lower_rate = as.vector(qmat_low),
           upper_rate = as.vector(qmat_high))
  
  # Format probabilities
  p_df <- as.data.frame(as.table(p_est)) |>
    rename(from = Var1, to = Var2, prob = Freq) |>
    mutate(lower_prob = as.vector(p_low),
           upper_prob = as.vector(p_high))
  
  # Combine
  left_join(q_df, p_df, by = c("from", "to"))
}


library(tictoc)
tic()
result_df <- prediction_grid |> 
  group_by(sex, age, int_date_decimal) |> 
  group_modify(~predict_hazard_and_prob(.x,
                                        .y,
                                        age_interval = .25,
                                        models = models)) |> 
  ungroup()
toc()

head(result_df)

result_df |> 
  mutate(from = from |> as.character() |> parse_number(),
         to = to |> as.character() |> parse_number(),
         transition = paste0("m",from,to)) |> 
  filter(to > from) |> 
  ggplot(aes(x=age,y=rate, color = transition, fill = transition)) +
  geom_line() +
  #geom_ribbon(mapping = aes(ymin = upper_rate, ymax = lower_rate), alpha = .3) +
  facet_wrap(int_date_decimal~sex) +
  scale_y_log10() +
  theme_minimal()
result_df |> filter(sex == "Male", from == "State 1", to == "State 2")
head(result_df)
