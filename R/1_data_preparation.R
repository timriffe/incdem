# ------------------------------------------------------------------- #
source("R/0_packages.R")
# ------------------------------------------------------------------- #
# Data
hrs <- read_csv("Data/riffe_incdem_20250522.csv")
# ------------------------------------------------------------------- #
# first pass processing
hrs_msm <- hrs |>
  # pick age range to fit to, based on plot of support.
  # if all ages 
  filter(between(age, 55, 97)) |> 
  mutate(
    int_date         = suppressWarnings(as.integer(int_date)),
    interview_date   = as_date(int_date, origin = "1960-01-01"),
    int_date_decimal = decimal_date(interview_date)
  ) |> 
  arrange(hhidpn, age) |>
  # mutate(age_diff = age - (int_date_decimal - decimal_date(as_date(birth_date, origin = "1960-01-01")))) |> 
  mutate(
    state_msm = state + 1  # msm expects states starting at 1
  ) |> 
  group_by(hhidpn) |> 
  # supposed to remove solitary observations.
  filter(n() > 1) |> 
  ungroup()
# ------------------------------------------------------------------- #
# all filtering before here; needed to create spline 
# basis with consistent age range in fit and in predictions.
# spline_basis_fit is also used later for this reason
spline_basis_fit   <- ns(hrs_msm$age, df = 3)
age_splines        <- as.data.frame(spline_basis_fit)
names(age_splines) <- paste0("age_spline", 1:3)
# ------------------------------------------------------------------- #
# This second part is supposed to eliminate impossible transitions.
hrs_msm <- hrs_msm |> 
  bind_cols(age_splines) |>
  group_by(hhidpn) |> 
  arrange(age) |> 
  mutate(
    ever_dementia = cumany(state == 1),
    state_clean   = case_when(
      state == 2                 ~ 2,  # preserve death
      ever_dementia & state == 0 ~ 1,  # impute "recovery" as still dementia
      TRUE ~ state
    ),
    died = cumany(state_clean == 2)
    ) |>
  filter(!(died & state_clean != 2)) |> 
  ungroup() |> 
  arrange(hhidpn, age) |> 
  # treat death times as exact, but othe transition
  # times as unknown.
  mutate(obstype = ifelse(state_msm == 3, 3, 1))
# ------------------------------------------------------------------- #
# RT: After cleaning we have to remove a solitary obs. one more time
# 190 persons, they provide no information
hrs_msm <- hrs_msm |>
  group_by(hhidpn) |> 
  # supposed to remove solitary observations.
  filter(n() > 1) |>
  ungroup()
# ------------------------------------------------------------------- #