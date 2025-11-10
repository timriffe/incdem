#!/usr/bin/env Rscript

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse
times <- suppressWarnings(as.integer(args[1]))
n_cores <- suppressWarnings(as.integer(args[2]))

# Validate
if (is.na(times) || is.na(n_cores)) {
  stop("Both arguments must be valid integers. Got: ",
       "times = ", args[1], ", n_cores = ", args[2])
}

if (times < 1 || n_cores < 1) {
  stop("Both times and n_cores must be positive integers. Got: ",
       "times = ", times, ", n_cores = ", n_cores)
}
if (length(times) != 1 || !is.numeric(times)) stop("`times` must be a single numeric value.")

source("R/00_package_and_functions.R")
source("R/01_prepare_hrs.R")
source("R/zzz_boot_test_functions.R")
source("R/zzz_boot_extra_functions.R")
options(group_boot.verbose = FALSE)

data <- read_csv("Data/hrs_to_fit_short.csv.gz", show_col_types = FALSE)  

# Define your Q matrix (adapt as needed)
Q <- rbind(
  c(-0.2, 0.1, 0.1),
  c(0, -.01, 0.1), 
  c(0, 0, 0)
)

# Run model
result <- fit_msm_boot(
  data = data,
  strat_vars = c("female"),
  covariate_var = c("period"),
  age_from_to = c(50, 100),
  age_int = 0.25,
  spline_df = NULL,
  spline_type = "ns",
  Q = Q,
  times = times,
  weight = "pwt",
  group = "hhidpn",
  n_cores = min(times, n_cores)
)

# Optionally write results
outfile <- glue("results/boot_result_t{times}_c{n_cores}.rds")
dir.create("results", showWarnings = FALSE)
saveRDS(result, outfile)
