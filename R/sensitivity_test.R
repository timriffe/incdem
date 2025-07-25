# ------------------------------------------------------------------- #
# load data
source("R/01_prepare_hrs.R")
# ------------------------------------------------------------------- #
# lets take a look at some results
# example 1: stratify by sex and use obs_date as a covariate
# NOTE: since covariate is continuous, I provide cont_grid
# that is used in the prediction part to predict data on these particular point
# in this example it will return predictions for 2000 2010 and 2020
# age grid is from 50 to 100
# natural splines with 3 df

results1 <- hrs_to_fit |> # .data argument for piping
  fit_msm(
    # main part
    strat_vars    = c("female"), # any vector of stratification variables 
    covariate_var = c("obs_date"), # any vector of covariates
    # age grid part
    age_from_to   = c(50, 100), # from first to second
    age_int       = 0.25,
    # continuous time grid
    cont_grid     = c(2000, 2010, 2020), # if covariate_var is not continuous set to NULL
    # splines part
    spline_df     = 3, # only for calc_spline = T
    spline_type   = "ns", # only for calc_spline = T
    calc_spline   = TRUE, # if false it does not matter what spline_df and spline_type are
    # CI part
    n_cores       = 4, # only used for ci_type = "bootstrap"
    # notice it is very low for testing
    B             = 2, # for ci_type = normal and ci_type = bootstrap
    ci_type       = "normal",
    conf_level    = 0.95, # also for normal and bootstrap
    # create Q matrix
    Q = rbind(
      c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, -.01,   0.1),  # dementia can go to death
      c(0, 0,   0)        # death is absorbing
    )
  )


# here is the corresponding plot for example 1
# We also have data on lower and upper CI so can plot them too
# NOTE: the CI will be very jagged since B = 2 is very low for testing speed
# set to something like 100 or better yet 1000 to make them smooth and nice
results1 |>
  # filter the rate option q (another option is probability - p)
  filter(type == "q") |>
  mutate(
    to = to |> as.character() |> parse_number(),
    from = from |> as.character() |> parse_number()
  ) |>
  filter(to > from) |>
  mutate(Period = as.factor(obs_date),
         transition = paste0("haz", from, to)) |>
  ggplot(aes(
    x        = age,
    y        = estimate, # here we can choose lower or upper too
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
  fit_msm(
    strat_vars    = c("female"),
    covariate_var = c("education", "stroke"),
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    cont_grid     = NULL, # !!!
    spline_df     = 3,
    spline_type   = "ns",
    calc_spline   = FALSE,
    n_cores       = 1,
    B             = 2,
    # create Q matrix
    Q = rbind(
      c(0, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, 0,   0.1),  # dementia can go to death
      c(0, 0,   0)     # death is absorbing
    )
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
  fit_msm(
    strat_vars    = c("female", "Period"),
    covariate_var = c("stroke"),
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    cont_grid     = NULL,
    spline_df     = 3,
    spline_type   = "ns",
    calc_spline   = FALSE,
    n_cores       = 1,
    B             = 2,
    # create Q matrix
    Q = rbind(
      c(0, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, 0,   0.1),  # dementia can go to death
      c(0, 0,   0)     # death is absorbing
    )
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
  fit_msm(
    strat_vars    = c("female"),
    covariate_var = c("stroke", "Period"),
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    cont_grid     = NULL,
    spline_df     = 3,
    spline_type   = "ns",
    calc_spline   = FALSE,
    n_cores       = 1,
    B             = 2,
    # create Q matrix
    Q = rbind(
      c(0, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, 0,   0.1),  # dementia can go to death
      c(0, 0,   0)     # death is absorbing
    )
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