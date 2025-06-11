

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

transitions_fit2 <-
  hrs_to_fit |>
  # before and after 2010 just as an example
  mutate(period = ifelse(obs_date < 2010, "<2010", "2011+")) |>
  # factor it
  mutate(period = as.factor(period)) |> 
  fit_msm_models(
    strat_vars    = c("female", "period"),
    covariate_var = NULL,
    age_int       = 0.25,
    cont_grid     = NULL,
    spline_df     = 3,
    spline_type   = "ns",
    age_from_to = c(50, 100),
    calc_spline   = TRUE,
    n_cores       = 1,
    B             = 2,
    # create Q matrix
    Q = rbind(
      c(0, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, 0,   0.1),  # dementia can go to death
      c(0, 0,   0)     # death is absorbing
    )
  )


transitions_fit2 |>
  mutate(to = to |> as.character() |> parse_number(),
         from = from |> as.character()|> parse_number(),
         transition = paste0("m",from,to),
         gender = female |> as.character(),
         gender = if_else(gender == "0","men","women")) |> 
  filter(to > from) |> 
  ggplot(aes(x        = age, 
             y        = rate, 
             color    = transition, 
             linetype = period)) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap( ~ gender) + 
  theme(legend.position = "bottom")


transitions_fit3 <-
  hrs_to_fit |>
  # before and after 2010 just as an example
  mutate(period = case_when(between(obs_date,2000,2006) ~ "period 1",
                            between(obs_date,2006,2012) ~ "period 2",
                            between(obs_date,2012,2020) ~ "period 3")) |>
  # factor it
  mutate(period = as.factor(period)) |> 
  fit_msm_models(
    strat_vars    = c("female", "period"),
    covariate_var = NULL,
    age_int       = 0.25,
    cont_grid     = NULL,
    spline_df     = 3,
    spline_type   = "ns",
    age_from_to = c(50, 100),
    calc_spline   = TRUE,
    n_cores       = 1,
    B             = 2,
    # create Q matrix
    Q = rbind(
      c(0, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, 0,   0.1),  # dementia can go to death
      c(0, 0,   0)     # death is absorbing
    )
  )



period_covariate3 <-
  transitions_fit3 |>
  mutate(to = to |> as.character() |> parse_number(),
         from = from |> as.character()|> parse_number()) |> 
  filter(to > from) |> 
  mutate(transition = paste0("haz", from, to)) |>
  ggplot(aes(x        = age, 
             y        = rate, 
             color    = transition, 
             linetype = period)) +
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

period_covariate3


