source("R/01_prepare_hrs.R")
source("R/zzz_boot_test_functions.R")
source("R/zzz_boot_extra_functions.R")
options(group_boot.verbose = FALSE)

# source("R/00_package_and_functions.R")
# hrs_to_fit_short <- hrs_to_fit |>
#   mutate(period = case_when(between(obs_date, 2000, 2006) ~ "period 1",
#                             between(obs_date, 2006, 2012) ~ "period 2",
#                             between(obs_date, 2012, 2020) ~ "period 3")) |>
#   select(hhidpn, female, education,
#          pwt = wtcrnh, age, hypertension:stroke,
#          state_msm, ever_dementia, obstype, period)
hrs_to_fit
hrs_to_fit |> 
  pull(obs_date)
hrs_to_fit_short2 <- hrs_to_fit |>
  mutate(period = case_when(between(obs_date, 2000, 2010) ~ "period 1",
                            between(obs_date, 2010, 2019) ~ "period 2",
                            TRUE ~ "other"),
         period = as.factor(period)) |>
  select(hhidpn, female, education,
         pwt = wtcrnh, age, hypertension:stroke,
         state_msm, ever_dementia, obstype, period) |> 
  filter(period != "other")

hrs_to_fit_short2 |> write_csv("Data/hrs_to_fit_short.csv.gz")

no_boot <- hrs_to_fit_short2 |>
  fit_msm(
    strat_vars    = c("female","period"),
    covariate_var = NULL,
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    spline_df     = 1,
    spline_type   = "ns",
    Q = rbind(
      c(-0.2, 0.1, 0.1),
      c(0,   -.01, 0.1),
      c(0,    0,   0)
    )
  )


no_boot |> 
  # mutate(from = from |> as.character() |> parse_number(),
  #        to = to |> as.character() |> parse_number()) |> 
  filter(to>from) |> 
  ggplot(aes(x=age,y=haz,color=interaction(from,to)))+
  geom_line() +
  facet_wrap(female~period)

library(tictoc)
tic()
booty <- hrs_to_fit_short2 |>
  fit_msm_boot(
    strat_vars    = c("female"),
    covariate_var = c("age","period"),
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    spline_df     = NULL,          
    spline_type   = "ns",
    Q = rbind(
      c(-0.2, 0.1, 0.1),
      c(0,   -.01, 0.1),
      c(0,    0,   0)
    ),
    times   = 4,
    weight  = "pwt",
    group   = "hhidpn",
    n_cores = 2
  )
toc()
booty$haz_u - booty$haz_l
booty |> 
  filter(to>from) |> 
  ggplot(aes(x=age,y=haz,color=interaction(from,to)))+
  geom_line()+
  geom_ribbon(aes(ymin=haz_l, ymax=haz_u, fill = interaction(from,to))) +
  facet_wrap(period~female)


library(rsample)

boot_rset <- group_bootstraps2(
  data   = hrs_to_fit_short2,
  group  = "hhidpn",
  weight = "pwt",
  times  = 1
)

sp <- boot_rset$splits[[1]]
d_resample <- rsample::analysis(sp)

# relabel duplicates exactly like the function does
d_resample <- d_resample %>%
  dplyr::group_by(hhidpn) %>%
  dplyr::mutate(.dup_id = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    hhidpn = dplyr::if_else(
      .dup_id == 1L,
      hhidpn,
      paste0(hhidpn, "_", .dup_id)
    )
  ) %>%
  dplyr::select(-.dup_id)

# Now test-fit
test_fit <- fit_msm(
  data         = d_resample,
  strat_vars   = c("female"),
  covariate_var= c("age","period"),
  age_from_to  = c(50,100),
  age_int      = 0.25,
  spline_df    = NULL,
  Q            = rbind(
    c(-0.2, 0.1, 0.1),
    c(0, -.01, 0.1),
    c(0, 0, 0)
  )
)

source("R/zzz_boot_extra_functions.R")  # to load msm_model()

# see what msm_model() is doing
debugonce(msm_model)

tmp_mod <- msm_model(
  data         = d_resample,
  Q            = rbind(
    c(-0.2, 0.1, 0.1),
    c(0, -.01, 0.1),
    c(0, 0, 0)
  ),
  covariate_var= c("age","period"),
  extra_covars = NULL
)
library(msm)

booty <- hrs_to_fit_short2 |>
  fit_msm_boot(
    strat_vars    = c("female"),
    covariate_var = c("age","period"),
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    spline_df     = NULL,
    spline_type   = "ns",
    Q = rbind(
      c(-0.2, 0.1, 0.1),
      c(0, -.01, 0.1),
      c(0, 0, 0)
    ),
    times   = 4,
    weight  = "pwt",
    group   = "hhidpn",
    n_cores = 2
  )
head(booty)
hrs_to_fit_short2 |> 
  filter(is.na(period))
earlier_attempt<-
fit_msm(hrs_to_fit_short2,
         strat_vars    = c("female","period"),
         covariate_var = NULL,
         age_from_to   = c(50, 100),
         age_int       = 0.25,
         spline_df     = 1,
         spline_type   = "ns",
         calc_spline = TRUE,
         Q = rbind(
           c(-0.2, 0.1, 0.1),
           c(0, -.01, 0.1),
           c(0, 0, 0)
         ),
         B   = 4,
        ci_type = "bootstrap",
         conf_level = 0.95,
         n_cores = 2)
earlier_attempt
earlier_attempt |> 
  mutate(from = parse_number(as.character(from)),
         to = parse_number(as.character(to))) |> 
  filter(to>from,
         type == "q") |> 
  ggplot(aes(x=age,y=estimate,color=interaction(from,to)))+
  geom_line()+
  # geom_ribbon(aes(ymin=lower, ymax=upper, fill = interaction(from,to))) +
  facet_wrap(period~female)
  

Q3 <- rbind(
  c(-0.2,  0.1, 0.1),
  c( 0.0, -0.01, 0.1),
  c( 0.0,  0.0, 0.0)
)

# A) build 3 stratified, group boot replicates (by female)
boot_tbl <- make_stratified_group_boot(
  data        = hrs_to_fit_short2,
  id_col      = "hhidpn",
  strata_vars = c("female"),
  times       = 3,
  seed        = 123
)

str(boot_tbl)
dim(boot_tbl)
inspect_bootstrap_reps(boot_tbl, id_col = "hhidpn", cat_covs = c("female","period"))
boot_tbl$data[[1]]

rep_results <- purrr::imap_dfr(boot_tbl$data[1:3], function(df, i) {
  cat("\nFitting replicate", i, "...\n")
  out <- fit_msm_on_dataframe(
    df,
    strat_vars    = c("female"),
    covariate_var = c("age","period"), # linear age + period effect, no splines
    age_from_to   = c(50,100),
    age_int       = 0.25,
    spline_df     = NULL,
    spline_type   = "ns",
    Q             = Q3
  )
  if (nrow(out)) out$replicate <- i
  out
})
rep_results %>%
  ungroup() |> 
  dplyr::filter(from == 1, to == 2) %>%
  dplyr::group_by(female, period, age) %>%
  dplyr::summarise(haz = mean(haz, na.rm = TRUE),
                   haz_l = min(haz, na.rm=TRUE),
                   haz_u = max(haz, na.rm = TRUE),.groups = "drop") 

rep_results |> 
  filter(female == 0,
         to>from) |> 
  ggplot(aes(x=age,y=haz,group = interaction(period, replicate), color = period)) +
  geom_line() +
  facet_wrap(to~from) +
  scale_y_log10()


