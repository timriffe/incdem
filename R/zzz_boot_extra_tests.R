source("R/01_prepare_hrs.R")
source("R/zzz_boot_test_functions.R")
source("R/zzz_boot_extra_functions.R")
options(group_boot.verbose = FALSE)

hrs_to_fit |> 
  mutate(date5 = obs_date - obs_date %% 5,
         year = floor(obs_date)) |> 
 pull(year) |> table()



hrs_to_fit_short2 <- hrs_to_fit |>
  mutate(period = case_when(between(obs_date, 2000, 2010) ~ "period 1",
                            between(obs_date, 2010, 2015) ~ "period 2",
                            obs_date > 2015 ~ "period 3"
                            TRUE ~ "other")) |>
  select(hhidpn, female, education,
         pwt = wtcrnh, age, hypertension:stroke,
         state_msm, ever_dementia, obstype, period) |> 
  filter(period != "other") |> 
  mutate(period = as.factor(period),
         year = floor(obs_date))

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


N <- 1000
n_cores <- 6
at_a_time <- n_cores * 3
loop_i <- ceiling(N / at_a_time)

for (i in 1:loop_i){
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
    times   = at_a_time,
    weight_col  = "pwt",
    id_col   = "hhidpn",
    n_cores = n_cores,
    parallel = "mclapply",
    return_replicates = TRUE
  )
toi   <- i * at_a_time
fromi <- toi - at_a_time + 1
addi  <- fromi - 1
booty$replicates$replicate <- booty$replicates$replicate + addi
namei <- paste0("replicates_",fromi,"-",toi,".csv.gz")
write_csv(booty$replicates, file = file.path("Data/boot_replicates", namei))
rm(booty); gc()
}

files <- dir("Data/boot_replicates")
booty <- vroom::vroom(file = file.path("Data/boot_replicates",files))
booty |> 
  filter(from < to) |> 
  group_by(period,female,from,to,age) |> 
  summarize(hazard = mean(haz),
            lower = quantile(haz,.025),
            upper = quantile(haz,.975)) |> 
  ggplot(aes(x=age,y=hazard, color = interaction(from, to))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = interaction(from, to)), alpha = .3, color = "transparent") +
  facet_wrap(period ~ female) +
  scale_y_log10()


test_future <- hrs_to_fit_short2 |>
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
    times   = 8,
    weight_col  = "pwt",
    id_col   = "hhidpn",
    n_cores = 4,
    parallel = "future",
    return_replicates = TRUE
  )
future::plan()
