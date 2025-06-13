source("R/00_package_and_functions.R")
source("R/01_prepare_hrs.R")

p3df4 <-
  hrs_to_fit |>
  # before and after 2010 just as an example
  mutate(period = case_when(between(obs_date,2000,2006) ~ "period 1",
                            between(obs_date,2006,2012) ~ "period 2",
                            between(obs_date,2012,2020) ~ "period 3")) |>
  # factor it
  mutate(period = as.factor(period)) |> 
  fit_msm(
    strat_vars    = c("female", "period"),
    age_int       = 0.25,
    spline_df     = 4,
    spline_type   = "ns",
    age_from_to = c(50, 100),
    calc_spline   = TRUE,
    n_cores       = 4,
    B             = 40,
    ci_type       = "normal",
    conf_level    = 0.95,
    # create Q matrix
    Q = rbind(
      c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, -.01,   0.1),  # dementia can go to death
      c(0, 0,   0)        # death absorbing
    )
  )
p3df3 <-
  hrs_to_fit |>
  # before and after 2010 just as an example
  mutate(period = case_when(between(obs_date,2000,2006) ~ "period 1",
                            between(obs_date,2006,2012) ~ "period 2",
                            between(obs_date,2012,2020) ~ "period 3")) |>
  # factor it
  mutate(period = as.factor(period)) |> 
  fit_msm(
    strat_vars    = c("female", "period"),
    age_int       = 0.25,
    spline_df     = 3,
    spline_type   = "ns",
    age_from_to = c(50, 100),
    calc_spline   = TRUE,
    n_cores       = 4,
    ci_type       = "normal",
    conf_level    = 0.95,
    # create Q matrix
    Q = rbind(
      c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, -.01,   0.1),  # dementia can go to death
      c(0, 0,   0)        # death absorbing
    )
  )


p3df2 <-
  hrs_to_fit |>
  # before and after 2010 just as an example
  mutate(period = case_when(between(obs_date,2000,2006) ~ "period 1",
                            between(obs_date,2006,2012) ~ "period 2",
                            between(obs_date,2012,2020) ~ "period 3")) |>
  # factor it
  mutate(period = as.factor(period)) |> 
  fit_msm(
    strat_vars    = c("female", "period"),
    age_int       = 0.25,
    spline_df     = 2,
    spline_type   = "ns",
    age_from_to = c(50, 100),
    calc_spline   = TRUE,
    n_cores       = 4,
    ci_type       = "normal",
    conf_level    = 0.95,
    # create Q matrix
    Q = rbind(
      c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, -.01,   0.1),  # dementia can go to death
      c(0, 0,   0)        # death absorbing
    )
  )


p3df1 <-
  hrs_to_fit |>
  # before and after 2010 just as an example
  mutate(period = case_when(between(obs_date,2000,2006) ~ "period 1",
                            between(obs_date,2006,2012) ~ "period 2",
                            between(obs_date,2012,2020) ~ "period 3")) |>
  # factor it
  mutate(period = as.factor(period)) |> 
  fit_msm(
    strat_vars    = c("female", "period"),
    age_int       = 0.25,
    spline_df     = 1,
    spline_type   = "ns",
    age_from_to = c(50, 100),
    calc_spline   = TRUE,
    n_cores       = 4,
    ci_type       = "normal",
    conf_level    = 0.95,
    # create Q matrix
    Q = rbind(
      c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, -.01,   0.1),  # dementia can go to death
      c(0, 0,   0)        # death absorbing
    )
  )

p3_linear <-
  hrs_to_fit |>
  # before and after 2010 just as an example
  mutate(period = case_when(between(obs_date,2000,2006) ~ "period 1",
                            between(obs_date,2006,2012) ~ "period 2",
                            between(obs_date,2012,2020) ~ "period 3")) |>
  # factor it
  mutate(period = as.factor(period)) |> 
  fit_msm(
    strat_vars    = c("female", "period"),
    age_int       = 0.25,
    age_from_to = c(50, 100),
    calc_spline   = FALSE,
    n_cores       = 4,
    ci_type       = "normal",
    conf_level    = 0.95,
    # create Q matrix
    Q = rbind(
      c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
      c(0, -.01,   0.1),  # dementia can go to death
      c(0, 0,   0)        # death absorbing
    )
  )
p3df4 <- 
  p3df4 |> 
  mutate(age_model = "df4")
p3df3 <- 
  p3df3 |> 
  mutate(age_model = "df3")
p3df2 <- 
  p3df2 |> 
  mutate(age_model = "df2")
p3df1 <- 
  p3df1 |> 
  mutate(age_model = "df1")
p3_linear <- 
  p3_linear |> 
  mutate(age_model = "linear")

compare_age <- bind_rows(p3df4,
                         p3df3,
                         p3df2,
                         p3df1,
                         p3_linear)

compare_age |> 
  filter(type == "q",
         from == "State 1",
         to == "State 2") |> 
  ggplot(aes(x=age,y=estimate, fill = age_model)) +
  geom_line(mapping=aes(color=age_model)) +
  facet_grid(vars(female),vars(period)) +
  theme_minimal() +
  scale_y_log10() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3)

p3df4 |> 
  filter(type == "q",
         from == "State 1",
         to == "State 2") |> 
  ggplot(aes(x=age,y=estimate)) +
  geom_line() +
  facet_grid(vars(female),vars(period)) +
  theme_minimal() +
  scale_y_log10() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3)

