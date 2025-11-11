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
hrs_to_fit_short2 <- hrs_to_fit |>
  mutate(period = case_when(between(obs_date, 2000, 2010) ~ "period 1",
                            between(obs_date, 2010, 2012) ~ "period 2")) |>
  select(hhidpn, female, education,
         pwt = wtcrnh, age, hypertension:stroke,
         state_msm, ever_dementia, obstype, period)
hrs_to_fit_short2 |> write_csv("Data/hrs_to_fit_short.csv.gz")

no_boot <- hrs_to_fit_short2 |>
  fit_msm(strat_vars    = c("female"), 
          covariate_var = c("period"), 
          age_from_to   = c(50, 100), 
          age_int       = 0.25,
          spline_df     = 2,    # only for calc_spline = T
          spline_type   = "ns", # only for calc_spline = T
          Q = rbind(
            c(-0.2, 0.1, 0.1),
            c(0, -.01,   0.1), 
            c(0, 0,   0) 
          ))


library(tictoc)
tic()
booty <- hrs_to_fit_short2 |>
  fit_msm_boot(
    strat_vars    = c("female"),
    covariate_var = c("period"),
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    spline_df     = 2,          # <--- turn splines ON
    spline_type   = "ns",
    Q = rbind(
      c(-0.2, 0.1, 0.1),
      c(0,   -.01, 0.1),
      c(0,    0,   0)
    ),
    times   = 10,
    weight  = "pwt",
    group   = "hhidpn",
    n_cores = 5
  )
toc()
no_boot$female |> unique()
booty |> 
  filter(to>from) |> 
  ggplot(aes(x=age,y=haz,color=interaction(from,to), linetype = female))+
  geom_line()+
  geom_ribbon(aes(ymin=haz_l, ymax=haz_u))
library(tibble)
library(purrr)
library(bench)
library(dplyr)

# Create a tibble of parameter combinations


# results |> 
#   mutate(n_cores = pmin(n_cores, times))

# peak_r_memory <- function(expr, log_interval = 1) {
  log_file <- tempfile(fileext = ".txt")
  logger_script <- "log_r_mem.sh"
  
  # Ensure logging script is in place
  if (!file.exists(logger_script)) {
    writeLines(c(
      "#!/bin/bash",
      "LOGFILE=$1",
      "> \"$LOGFILE\"",
      "while true; do",
      "  ts=$(date +%s)",
      "  mem=$(ps -e -o comm,rss | grep -E 'R$|Rscript$|RStudio' | awk '{sum += $2} END {print sum}')",
      "  echo \"$ts $mem\" >> \"$LOGFILE\"",
      paste0("  sleep ", log_interval),
      "done"
    ), logger_script)
    Sys.chmod(logger_script, mode = "0755")
  }
  
  # Start memory logging in background
  pid <- system(glue::glue("nohup ./{logger_script} {log_file} & echo $!"), intern = TRUE)
  logger_pid <- as.integer(pid)
  
  on.exit({
    # Ensure the logger is stopped even on error
    system(glue::glue("kill {logger_pid}"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  }, add = TRUE)
  
  # Run target expression and time it
  result <- eval(expr)
  
  # Stop the logger
  system(glue::glue("kill {logger_pid}"))
  
  # Read memory log
  mem_log <- readr::read_table(log_file, col_names = c("timestamp", "rss_kb"), show_col_types = FALSE)
  peak_mem_mb <- max(mem_log$rss_kb, na.rm = TRUE) / 1024
  
  list(result = result, peak_memory_mb = peak_mem_mb)
}
# 
# times <- c(2, 4, 8, 16)
# n_cores <- c(1, 2, 4)
# grid <- crossing(
#   times = times,
#   n_cores = n_cores
# ) |> 
#   filter(times >= n_cores)
# 
# # Initialize results tibble
# results <- grid %>%
#   mutate(
#     metrics = pmap(list(times, n_cores), function(times, n_cores) {
#       # Measure time and memory
#       start_time <- Sys.time()
#       run <- peak_r_memory({
#         fit_msm_boot(
#           data          = hrs_to_fit_short2,
#           strat_vars    = c("female"), 
#           covariate_var = c("period"), 
#           age_from_to   = c(50, 100), 
#           age_int       = 0.25,
#           spline_df     = NULL,
#           spline_type   = "ns",
#           Q             = rbind(
#             c(-0.2, 0.1, 0.1),
#             c(0, -.01, 0.1), 
#             c(0, 0, 0)
#           ),
#           times         = times,
#           weight        = "pwt",
#           group         = "hhidpn",
#           n_cores       = min(times, n_cores)
#         )
#       })
#       end_time <- Sys.time()
#       
#       tibble::tibble(
#         median_time = end_time - start_time,
#         max_mem     = as.numeric(run$peak_memory_mb, units = "MB")
#       )
#     })
#   ) %>%
#   tidyr::unnest(metrics)

