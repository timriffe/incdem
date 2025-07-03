source("R/01_prepare_hrs.R")
source("R/zzz_boot_test_functions.R")

hrs_to_fit_short <- hrs_to_fit |>
  mutate(period = case_when(between(obs_date, 2000, 2006) ~ "period 1",
                            between(obs_date, 2006, 2012) ~ "period 2",
                            between(obs_date, 2012, 2020) ~ "period 3"),
         # important
         period = as.factor(period)) |>
  select(hhidpn, female, education,
         wtcrnh, age, hypertension:stroke,
         state_msm, ever_dementia, obstype, period)

# two examples with and with no boot
no_boot <- hrs_to_fit_short |>
  fit_msm(strat_vars    = c("female"), # any vector of stratification variables 
          covariate_var = c("period"), # any vector of covariates
          # age grid part
          age_from_to   = c(50, 100), # from first to second
          age_int       = 0.25,
          # continuous time grid
          cont_grid     = NULL, # if covariate_var is not continuous set to NULL
          # splines part
          spline_df     = 3, # only for calc_spline = T
          spline_type   = "ns", # only for calc_spline = T
          calc_spline   = TRUE, # if false it does not matter what spline_df and spline_type are
          # CI part
          n_cores       = 8, # only used if ci = TRUE
          # notice it is very low for testing
          B             = 5, # number of bootstraps if FALSE it is ignored
          ci            = FALSE, # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          conf_level    = 0.95, # also for bootstrap
          # create Q matrix
          Q = rbind(
            c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
            c(0, -.01,   0.1),  # dementia can go to death
            c(0, 0,   0)        # death is absorbing
          ),
          # weather prevalence calculation is needed
          calc_prev = TRUE)

# now with CI
with_boot <- hrs_to_fit_short |>
  fit_msm(strat_vars    = c("female"), # any vector of stratification variables 
          covariate_var = c("period"), # any vector of covariates
          # age grid part
          age_from_to   = c(50, 100), # from first to second
          age_int       = 0.25,
          # continuous time grid
          cont_grid     = NULL, # if covariate_var is not continuous set to NULL
          # splines part
          spline_df     = 3, # only for calc_spline = T
          spline_type   = "ns", # only for calc_spline = T
          calc_spline   = TRUE, # if false it does not matter what spline_df and spline_type are
          # CI part
          n_cores       = 8, # only used for ci_type = "bootstrap"
          # notice it is very low for testing
          B             = 5, # for ci_type = normal and ci_type = bootstrap
          ci            = TRUE, # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          conf_level    = 0.95, # also for normal and bootstrap
          # create Q matrix
          Q = rbind(
            c(-0.2, 0.1, 0.1),  # healthy can go to dementia or death
            c(0, -.01,   0.1),  # dementia can go to death
            c(0, 0,   0)        # death is absorbing
          ),
          # weather prevalence calculation is needed
          calc_prev = TRUE)


# ---------------------------------------------------------------------------- #
# here is the mortality adjustment for estimate of mu
# NOTE Prevalence now is only calculates for the estimate and not for 
# CI High and Low. I can add that easily though if we approve the routine
# just mutate across
# HMD data for ages 50:100 is essentially a line (Gompertz)
# NOTE: This file is created from Dx and Ex from HMD with script extra step HMD
load("Data/hmd_mort.RData")

hmd_25 <- hmd_25 |> 
  select(-qx) |>
  mutate(female = as.factor(female))
# ---------------------------------------------------------------------------- #
# now that we have prevalence mx and hmd data
# we calculate Ra and apply formula 5 to point estimates
result <- with_boot |>
  select(female, period, age, from, to, estimate, prevalence) |>
  filter(to == "State 3") |>
  mutate(to = "D") |>
  mutate(from = ifelse(from == "State 1", "H", "U")) |>
  unite("trns", c(from, to), sep = "-") |>
  pivot_wider(names_from  = trns,
              values_from = c(estimate, prevalence)) |>
  select(-`prevalence_H-D`, prevalence = `prevalence_U-D`,
         `H-D` = `estimate_H-D`, `U-D` = `estimate_U-D`) |>
  full_join(hmd_25) |>
  # rate ratio
  mutate(Ra = `U-D` / `H-D`) |>
  group_by(female, period) |>
  # new rates
  mutate(
    mh_new = mu / (1 - prevalence + prevalence * Ra), # formula 5
    mu_new = mh_new * Ra) |>
  # just wrangling
  select(female, period, age, 
         mh_old = `H-D`, 
         mu_old = `U-D`,
         mh_new,
         mu_new) |>
  ungroup() |>
  pivot_longer(cols = c(mu_old, mu_new, mh_old, mh_new),
               names_to = "type", values_to = "value") |>
  separate(type, into = c("rate", "version"), sep = "_")

# there is a difference expecially in mu
# males
result |>
  filter(female == 0) |>
  ggplot(aes(x = age, y = value, color = version)) +
  geom_line() +
  facet_grid(rate ~ period) +
  labs(title = "Old vs New Mortality Rates",
       y = "Rate",
       color = "Version") +
  theme_minimal() + 
  scale_x_log10() + 
  theme_bw() + 
  theme(legend.position = "bottom")

# and females
result |>
  filter(female == 1) |>
  ggplot(aes(x = age, y = value, color = version)) +
  geom_line() +
  facet_grid(rate ~ period) +
  labs(title = "Old vs New Mortality Rates",
       y = "Rate",
       color = "Version") +
  theme_minimal() + 
  scale_x_log10() + 
  theme_bw() + 
  theme(legend.position = "bottom")