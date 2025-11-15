
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

# ----------------------------------------------------
# (1) age linear, 5-year period as strata
# ----------------------------------------------------
N <- 1000
n_cores <- 6
at_a_time <- n_cores * 3
loop_i <- ceiling(N / at_a_time)

for (i in 1:loop_i){
  booty <- hrs_to_fit_prepped |>
    fit_msm_boot(
      strat_vars    = c("female","period5"),
      covariate_var = c("age"),
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
  write_csv(booty$replicates, file = file.path("Data/boot_replicates1", namei))
  rm(booty); gc()
}

# ----------------------------------------------------
# (2) age linear, year linear
# ----------------------------------------------------
N <- 1000
n_cores <- 6
at_a_time <- n_cores * 3
loop_i <- ceiling(N / at_a_time)

for (i in 1:loop_i){
  booty <- hrs_to_fit_prepped |>
    fit_msm_boot(
      strat_vars    = c("female"),
      covariate_var = c("age","year"),
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
  write_csv(booty$replicates, file = file.path("Data/boot_replicates2", namei))
  rm(booty); gc()
}


# ----------------------------------------------------
# (3) age spline df2, year linear
# ----------------------------------------------------
N <- 1000
n_cores <- 6
at_a_time <- n_cores * 2
loop_i <- ceiling(N / at_a_time)

for (i in 1:loop_i){
  booty <- hrs_to_fit_prepped |>
    fit_msm_boot(
      strat_vars    = c("female"),
      covariate_var = c("year"),
      age_from_to   = c(50, 100),
      age_int       = 0.25,
      spline_df     = 2,
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
  write_csv(booty$replicates, file = file.path("Data/boot_replicates3", namei))
  rm(booty); gc()
}


# files <- dir("Data/boot_replicates1")
# booty <- vroom::vroom(file = file.path("Data/boot_replicates1",files))