
# --------------------
# Here we should estimate the empirical prevalence,
# as opposed to the stationary prevalence. This still needs
# to be modeled, and bootstrapped to combine w the initial
# hazards estimate at the adjustment step.
# --------------------


source("R/00_dependencies.R")
source("R/01_functions.R")
source("R/02_prepare_hrs.R")
source("R/01_functions.R")
hrs_to_fit_prepped <- 
  hrs_to_fit |>
  mutate(period5 = case_when(between(obs_date, 2004, 2010) ~ "period 1",
                             between(obs_date, 2010, 2015) ~ "period 2",
                             obs_date > 2015 ~ "period 3",
                             TRUE ~ "other"),
         year = floor(obs_date)) |>
  select(hhidpn, female, education,
         pwt = wtcrnh, age, hypertension:stroke,
         state_msm, ever_dementia, obstype, period5,
         year) |> 
  filter(period5 != "other",
         year > 2003,
         state_msm != 3) |> 
  mutate(period = as.factor(period5))

set.seed(1981)
boot_design <- group_bootstraps2(
  data   = hrs_to_fit_prepped,
  group  = "hhidpn",
  times  = 1000,
  weight = "pwt"
)
# -----------------------
# model 1 prevalence
# -----------------------
# increasing prevalence!
booty1 <-
  hrs_to_fit_prepped |> 
  fit_prev_boot(
  strat_vars    = c("female","period5"),
  covariate_var = "age",
  age_from_to   = c(50, 100),
  age_int       = 0.25,
  condition_state = 2L,
  spline_df     = NULL,
  spline_type   = "ns",
  state_var = "state_msm",
  weight_col  = "pwt",
  id_col   = "hhidpn",
  times           = 1008,
  ci_level        = 0.95,
  return_replicates = TRUE,
  boot_rset = boot_design
) %>% '[['("replicates")
write_csv(booty1, file = "Data/model1/prev_replicates.csv.gz")

 booty1 |> 
   group_by(female, period5, age) |> 
   summarize(prev_median = median(prevalence),
             lower= quantile(prevalence, .025),
             upper = quantile(prevalence, .975), .groups="drop") |> 
   ggplot(aes(x = age, y = prev_median, color = period5)) +
   geom_line() +
   geom_ribbon(aes(ymin=lower, ymax=upper, fill=period5),alpha=.3,color="transparent")+
   facet_wrap(~female)
  
# -----------------------
# model 2 prevalence
# -----------------------
# decreasing prevalence!
booty2 <- 
  hrs_to_fit_prepped |> 
  fit_prev_boot(
  strat_vars    = c("female"),
  covariate_var = c("age","year"),
  age_from_to   = c(50, 100),
  age_int       = 0.25,
  condition_state = 2L,
  spline_df     = NULL,
  spline_type   = "ns",
  state_var = "state_msm",
  weight_col  = "pwt",
  id_col   = "hhidpn",
  times           = 1008,
  ci_level        = 0.95,
  return_replicates = TRUE,
  boot_rset = boot_design
)  %>% '[['("replicates")
 write_csv(booty2, file = "Data/model2/prev_replicates.csv.gz")
 
 booty2 |> 
   mutate(year = as.integer(year)) |> 
   filter(year %% 5 == 0) |> 
   group_by(female, year, age) |> 
   summarize(prev_median = median(prevalence),
             lower= quantile(prevalence, .025),
             upper = quantile(prevalence, .975), .groups="drop") |> 
   ggplot(aes(x = age, y = prev_median, color = as.factor(year))) +
   geom_line() +
   geom_ribbon(mapping=aes(ymin=lower, ymax=upper, fill=as.factor(year)),
               alpha=.3,color="transparent")+
   facet_wrap(~female)

 booty3 <- 
   hrs_to_fit_prepped |> 
   fit_prev_boot(
     strat_vars    = c("female"),
     covariate_var = c("age","year"),
     age_from_to   = c(50, 100),
     age_int       = 0.25,
     condition_state = 2L,
     spline_df     = NULL,
     spline_type   = "ns",
     state_var = "state_msm",
     weight_col  = "pwt",
     id_col   = "hhidpn",
     times           = 1008,
     ci_level        = 0.95,
     return_replicates = TRUE,
     boot_rset = boot_design
   )  %>% '[['("replicates")
 write_csv(booty3, file = "Data/model3/prev_replicates.csv.gz")
rm(booty1);rm(booty2);rm(booty3);gc()









# Furter possible specifications:
do_this <- FALSE
if (do_this){
 # point estimates ok...
 hrs_to_fit_prepped |> 
fit_prev(
  strat_vars    = c("female"),
  covariate_var = c("year"),
  age_from_to   = c(50, 100),
  age_int       = 0.25,
  condition_state = 2L,
  spline_df     = 2,
  spline_type   = "ns",
  state_var = "state_msm",
  weight_col  = "pwt",
  id_col   = "hhidpn"
)|> 
  filter(year %% 5 == 0) |> 
  ggplot(aes(x = age, y = prevalence, color = as.factor(year))) +
  geom_line() +
  facet_wrap(~female)


hrs_to_fit_prepped |> 
  fit_prev(
    strat_vars    = c("female"),
    covariate_var = c("year"),
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    condition_state = 2L,
    spline_df     = 3,
    spline_type   = "ns",
    state_var = "state_msm",
    weight_col  = "pwt",
    id_col   = "hhidpn"
  )|> 
  filter(year %% 5 == 0) |> 
  ggplot(aes(x = age, y = prevalence, color = as.factor(year))) +
  geom_line() +
  facet_wrap(~female)


hrs_to_fit_prepped |> 
  fit_prev(
    strat_vars    = c("female"),
    covariate_var = c("year"),
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    condition_state = 2L,
    spline_df     = 3,
    spline_type   = "ns",
    state_var = "state_msm",
    weight_col  = "pwt",
    id_col   = "hhidpn"
  )|> 
  filter(year %% 5 == 0) |> 
  ggplot(aes(x = age, y = prevalence, color = as.factor(year))) +
  geom_line() +
  facet_wrap(~female)
}