
# --------------------
# Here we should estimate the empirical prevalence,
# as opposed to the stationary prevalence. This still needs
# to be modeled, and bootstrapped to combine w the initial
# hazards estimate at the adjustment step.
# --------------------


source("R/00_package_and_functions.R")
source("R/01_prepare_hrs.R")
source("R/00_package_and_functions.R")
hrs_to_fit_prepped <- 
  hrs_to_fit |>
  mutate(period5 = case_when(between(obs_date, 2000, 2010) ~ "period 1",
                             between(obs_date, 2010, 2015) ~ "period 2",
                             obs_date > 2015 ~ "period 3",
                             TRUE ~ "other"),
         year = floor(obs_date)) |>
  select(hhidpn, female, education,
         pwt = wtcrnh, age, hypertension:stroke,
         state_msm, ever_dementia, obstype, period5,
         year) |> 
  filter(period5 != "other",
         year > 2003) |> 
  mutate(period = as.factor(period5))

# increasing prevalence!
hrs_to_fit_prepped |> 
fit_prev(
  strat_vars    = c("female","period5"),
  covariate_var = c("age"),
  age_from_to   = c(50, 100),
  age_int       = 0.25,
  condition_state = 2L,
  spline_df     = NULL,
  spline_type   = "ns",
  state_var = "state_msm",
  #weight_col  = "pwt",
  id_col   = "hhidpn"
) |> 
  ggplot(aes(x = age, y = prevalence, color = period5)) +
  geom_line() +
  facet_wrap(~female)
  






