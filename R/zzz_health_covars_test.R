source("R/01read_and_prepare_hrs.R")
# source("R/01a_read_recode_hrs.R")
# source("R/01b_data_preparation.R")




# at this time, health variables are all time-varying ever-had variables

# hrs_to_fit |> 
#   group_by(hhidpn) |> 
#   summarize(n_st = unique(stroke) |> length()) |> 
#   pull(n_st) |> unique()


# (1) "hypertension"
fitted_hyp_strata <-
  hrs_to_fit |> 
  fit_msm(strat_vars    = c("female", "hypertension"),
                      age_int       = 0.25,
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
                 )) |> 
  mutate(health_var = "hypertension",
         type = "strata")

fitted_hyp_cov <-
  hrs_to_fit |> 
  fit_msm(strat_vars    = c("female"),
                      covariate_var = c("hypertension"),
                      age_int       = 0.25,
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
                 )) |> 
  mutate(health_var = "hypertension",
         type = "covariate")

# (2) "diabetes"
  
  fitted_diabetes_strata <-
    hrs_to_fit |> 
    fit_msm(strat_vars    = c("female", "diabetes"),
                        age_int       = 0.25,
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
                   )) |> 
    mutate(health_var = "diabetes",
           type = "strata")
  
  fitted_diabetes_cov <-
    hrs_to_fit |> 
    fit_msm(strat_vars    = c("female"),
                        covariate_var = c("diabetes"),
                        age_int       = 0.25,
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
                   )) |> 
    mutate(health_var = "diabetes",
           type = "covariate")
  
# (3) "heart_disease"
  
  fitted_heart_disease_strata <-
    hrs_to_fit |> 
    fit_msm(strat_vars    = c("female", "heart_disease"),
                        age_int       = 0.25,
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
                   )) |> 
    mutate(health_var = "heart_disease",
           type = "strata")
  
  fitted_heart_disease_cov <-
    hrs_to_fit |> 
    fit_msm(strat_vars    = c("female"),
                        covariate_var = c("heart_disease"),
                        age_int       = 0.25,
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
                   )) |> 
    mutate(health_var = "heart_disease",
           type = "covariate")
  
# (4) "stroke"
  fitted_stroke_strata <-
    hrs_to_fit |> 
    fit_msm(strat_vars    = c("female", "stroke"),
                        age_int       = 0.25,
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
                   )) |> 
    mutate(health_var = "stroke",
           type = "strata")
  
  fitted_stroke_cov <-
    hrs_to_fit |> 
    fit_msm(strat_vars    = c("female"),
                        covariate_var = c("stroke"),
                        age_int       = 0.25,
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
                   )) |> 
    mutate(health_var = "stroke",
           type = "covariate")
  
# bind all results:
fitted_hyp <-
  fitted_hyp_strata |> 
  bind_rows(fitted_hyp_cov) |> 
  rename(status = hypertension)
fitted_diabetes <-
  fitted_diabetes_strata |> 
  bind_rows(fitted_diabetes_cov) |> 
  rename(status = diabetes)
fitted_heart_disease <-
  fitted_heart_disease_strata |> 
  bind_rows(fitted_heart_disease_cov) |> 
  rename(status = heart_disease)
fitted_stroke <-
  fitted_stroke_strata |> 
  bind_rows(fitted_stroke_cov) |> 
  rename(status = stroke)
  
health_variable_tests <- 
bind_rows(fitted_hyp, 
          fitted_diabetes, 
          fitted_heart_disease, 
          fitted_stroke) |> 
  mutate(to = to |> as.character() |> parse_number(),
         from = from |> as.character()|> parse_number(),
         transition = paste0("m",from,to),
         gender = female |> as.character(),
         gender = if_else(gender == "0","men","women")) |> 
  select(-female)

write_csv(health_variable_tests, "Data/health_var_tests.csv.gz")


health_variable_tests |> 
  filter(transition == "m12") |> 
  ggplot(aes(x = age, y = rate, color = type, linetype = status)) +
  geom_line() +
  scale_y_log10() +
  facet_grid(vars(gender),vars(health_var)) +
  theme_minimal()





