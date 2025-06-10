# ------------------------------------------------------------------- #
# load data
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


# ------------------------------------------------------------------- #
# create obs_date code by TR
hrs_to_fit <- hrs_msm |>
  # a test condition to check pandemic effect on mort trend
  mutate(
    birth_date         = as_date(birth_date, origin = "1960-01-01"),
    birth_date_decimal = decimal_date(birth_date),
    obs_date           = birth_date_decimal + age
  ) |>
  filter(obs_date < (as_date("2019-dec-31") |> decimal_date())) |>
  # factor variables
  mutate(across(c(female, race, hispanic, education,
                  hypertension, diabetes, heart_disease,
                  stroke, ever_dementia), ~ as.factor(.)))

# ------------------------------------------------------------------- #
# lets take a look at some results
# example 1: stratify by sex and use obs_date as a covariate
# NOTE: since covariate is continuous, I provide cont_grid
# that is used in the prediction part to predict data on these particular point
# in this example it will return predictions for 2000 2010 and 2020
# age grid is from 50 to 100
# natural splines with 3 df
results1 <- hrs_to_fit |>
  fit_msm_sensitivity(
    # main part
    strat_vars    = c("female"),
    covariate_var = c("obs_date"),
    # age grid part
    age_pred_grid = c(50, 100),
    age_int       = 0.25,
    # continuous time grid
    cont_grid     = c(2000, 2010, 2020),
    # splines part
    spline_df     = 3,
    spline_type   = "ns"
  )

# here is the corresponding plot for example 1
results1 |>
  mutate(
    to = to |> as.character() |> parse_number(),
    from = from |> as.character() |> parse_number()
  ) |>
  filter(to > from) |>
  mutate(Period = as.factor(obs_date),
         transition = paste0("haz", from, to)) |>
  ggplot(aes(
    x        = age,
    y        = rate,
    color    = transition,
    linetype = Period
  )) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap( ~ female) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------- #
# Example 2
# stratify by sex and use stroke and education as covariate
# NOTE: since there is no continuous covariates  cont_grid is set to NULL
results2 <- hrs_to_fit |>
  # filter(!is.na(Date), !is.na(hypertension)) |>
  fit_msm_sensitivity(
    strat_vars    = c("female"),
    covariate_var = c("education", "stroke"),
    age_pred_grid = c(50, 100),
    age_int       = 0.25,
    cont_grid     = NULL, # !!!
    spline_df     = 3,
    spline_type   = "ns"
  )

# here is the figure for example 2
results2 |>
  mutate(
    to = to |> as.character() |> parse_number(),
    from = from |> as.character() |> parse_number()
  ) |>
  filter(to > from) |>
  mutate(transition = paste0("haz", from, to)) |>
  ggplot(aes(
    x        = age,
    y        = rate,
    color    = transition,
    linetype = stroke
  )) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap(education ~ female) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------- #
# Example 3
# Lets separate data by time in 2 periods
results3 <- hrs_to_fit |>
  # before and after 2010 just as an example
  mutate(Period = ifelse(obs_date < 2010, "<2010", "2011+")) |>
  # factor it
  mutate(Period = as.factor(Period)) |> 
# Here we use stroke as covariate and sex and Period as stratas
  # filter(!is.na(Date), !is.na(hypertension)) |>
  fit_msm_sensitivity(
    strat_vars    = c("female", "Period"),
    covariate_var = c("stroke"),
    age_pred_grid = c(50, 100),
    age_int       = 0.25,
    cont_grid     = NULL,
    spline_df     = 3,
    spline_type   = "ns"
  )

# here is the plot for example 3
results3 |>
  mutate(to = to |> as.character() |> parse_number(),
         from = from |> as.character()|> parse_number()) |> 
  filter(to > from) |> 
  mutate(transition = paste0("haz", from, to)) |>
  ggplot(aes(x        = age, 
             y        = rate, 
             color    = transition, 
             linetype = Period)) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap(stroke ~ female) + 
  theme(legend.position = "bottom")

# ------------------------------------------------------------------- #
# Example 4
# same as 3 but now Period is in covariates
results4 <- hrs_to_fit |>
  # filter(!is.na(Date), !is.na(hypertension)) |>
  fit_msm_sensitivity(
    strat_vars    = c("female"),
    covariate_var = c("stroke", "Period"),
    age_pred_grid = c(50, 100),
    age_int       = 0.25,
    cont_grid     = NULL,
    spline_df     = 3,
    spline_type   = "ns"
  )

# here is the plot for example 4
results4 |>
  mutate(to = to |> as.character() |> parse_number(),
         from = from |> as.character()|> parse_number()) |> 
  filter(to > from) |> 
  mutate(transition = paste0("haz", from, to)) |>
  ggplot(aes(x        = age, 
             y        = rate, 
             color    = transition, 
             linetype = Period)) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap(stroke ~ female) + 
  theme(legend.position = "bottom")

# ------------------------------------------------------------------- #
# end