# Load required packages
source("R/00_package_and_functions.R")

# Select necessary variables from main dataset
max_wave <- 16
hrs_file <- if_else(max_wave == 16, "randhrs1992_2022v1.dta","randhrs1992_2020v2.dta")
vars <- c("hhidpn", "hhid", "pn", "hacohort", "rabdate", "raddate",
          paste0("r", 5:max_wave, "iwstat"),
          paste0("r", 5:max_wave, "iwmid"),
          paste0("r", 5:max_wave, "hibpe"),
          paste0("r", 5:max_wave, "diabe"),
          paste0("r", 5:max_wave, "hearte"),
          paste0("r", 5:max_wave, "stroke"),
          "ragender", "raracem", "rahispan", "rabyear", "raeduc")

# Load main RAND HRS dataset
hrs_in <- read_dta(file.path("Data",hrs_file), 
                   encoding = "UTF-8", 
                   col_select = all_of(vars)) |> 
  zap_labels() |> 
  mutate(
    hhid = as.character(hhid),
    pn = as.character(pn),
    hhidpn = paste0(hhid, pn)
  )

# Load cognition dataset
cog_colnames <- names(read_dta("Data/cogfinalimp_9520wide.dta", n_max = 0))
cog_cols_needed <- c("hhid", "pn", grep("^cogfunction\\d{4}$", cog_colnames, value = TRUE))

cog <- read_dta("./Data/cogfinalimp_9520wide.dta", 
                encoding = "UTF-8", 
                col_select = all_of(cog_cols_needed)) |>
  zap_labels() |>
  mutate(
    hhid = as.character(hhid),
    pn = as.character(pn),
    hhidpn = paste0(hhid, pn)
  ) |>
  select(hhidpn, starts_with("cogfunction"))

# Pre-process HRS dataset
hrs_long <- 
  hrs_in |>
  mutate(
    female = as.integer(ragender == 2),
    race = case_when(
      raracem == 1 ~ "white",
      raracem == 2 ~ "black",
      raracem == 3 ~ "other"
    ),
    hispanic = if_else(rahispan == 1, "hispanic", "non-hispanic"),
    education = case_when(
      raeduc %in% c(1, 2) ~ "ls hs/ged",
      raeduc == 3 ~ "hs",
      raeduc %in% c(4, 5) ~ ">hs"
    ),
    birth_year = rabyear,
    cohort = hacohort
  ) |>
  pivot_longer(
    cols = matches("^r\\d{1,2}(iwstat|iwmid|hibpe|diabe|hearte|stroke)$"),
    names_to = c("wave", ".value"),
    names_pattern = "r(\\d{1,2})([a-z]+)"
  ) |>
  mutate(wave = as.integer(wave)) |>
  filter(iwstat %in% c(1, 4, 5)) |>
  mutate(
    age_interview = (iwmid - rabdate) / 365.25,
    age_death     = (raddate - rabdate) / 365.25,
    age           = if_else(iwstat == 5, age_death, age_interview),
    int_date      = if_else(iwstat == 5, raddate, iwmid)
  ) |>
  group_by(hhidpn) |> 
  mutate(n=n()) |> 
  ungroup() |> 
  filter(n > 1) |> 
  arrange(hhidpn, wave) |>
  group_by(hhidpn) |>
  mutate(age_imputed = impute_age(age, wave)) |> 
  ungroup() |>
  group_by(hhidpn) |>
  mutate(across(c(hibpe, diabe, hearte, stroke),
                ~ cummax(replace_na(., 0)),
                .names = "{.col}_status")) |>
  ungroup()

# stata cage# stata code:
# * recode 1/2 = 0 (no dementia), 3 = 1 (dementia)
# recode cogfunction (1/2 = 0) (3 = 1), gen(cog_dementia)
# 

# Cognition: reshape and process
cog_long <- 
  cog |>
  pivot_longer(
    cols = starts_with("cogfunction"),
    names_to = "year",
    names_pattern = "cogfunction(\\d{4})",
    values_to = "cogfunction"
  ) |>
  mutate(year = as.integer(year)) |>
  # ignore cog status before 2000; inconsistent measurement
  # until 2000+
  # but note, there are 2335 cases of a 'first' dementia observed
  # before 2000. a possible sensitivity analyses
  filter(year >= 2000) |> 
  arrange(hhidpn,year) |> 
  # Recode: 1/2 = 0, 3 = 1
  mutate(cog_dementia = case_when(
    cogfunction %in% 1:2 ~ 0L,
    cogfunction == 3     ~ 1L,
    TRUE                 ~ NA_integer_)) |> 
  arrange(hhidpn, year) |>
  group_by(hhidpn) |> 
  # downup rather than updown: has to do with row ordering not time.
  tidyr::fill(cog_dementia,.direction = "downup") |> 
  mutate(state = cummax(cog_dementia)) |> 
  
  # Forward fill
  # arrange(hhidpn, year) |>
  # group_by(hhidpn) |>
  # mutate(cog_dementia = zoo::na.locf(cog_dementia, na.rm = FALSE)) |>
  # 
  # # Backward fill
  # arrange(hhidpn, desc(year)) |>
  # mutate(cog_dementia = zoo::na.locf(cog_dementia, fromLast = TRUE, na.rm = FALSE)) |>
  # arrange(hhidpn, year) |>
  #   # mutate(
  #   #dementia = cummax(replace_na(cog_dementia, 0)),
  #   state = dementia
  # ) |>
  # Running max and state

  ungroup()

# Map year to wave (2000 = wave 5, ..., 2020 = wave 15)
year_wave_map <- tibble(
  year = seq(2000, 2022, by = 2),
  wave = 5:16
)

cog_long <- cog_long |>
  left_join(year_wave_map, by = join_by(year)) 

# Merge cognition into HRS data and assign final state
hrs_joined <- hrs_long |>
  left_join(cog_long, by = c("hhidpn", "wave")) |>
  mutate(
    state = case_when(
      iwstat == 5              ~ 2L,         # deceased
      !is.na(state)            ~ state,      # dementia status
      TRUE                     ~ NA_integer_ # NA otherwise
    )
  ) |> 
  select(-stroke, -hibpe, -diabe, -hearte, -cogfunction, -cog_dementia, -n) |> 
  rename(birth_date = rabdate,
         death_date = raddate,
         hypertension = hibpe_status,
         diabetes = diabe_status,
         heart_disease = hearte_status,
         stroke = stroke_status) |> 
  mutate(age = age_imputed) |> 
  select(-age_imputed)

# View result
# write_csv(hrs_joined, "./Data/rand_hrs_processed.csv")
 # hrs_jordan <- read_csv("Data/riffe_incdem_20250522.csv") |> 
 #   mutate(hhidpn = sprintf("%09.0f", hhidpn)) |> 
 #   select(hhidpn, wave, age_jordan = age, state_jordan = state)
 # 
 # hrs_joined |> 
 #   select(hhidpn,
 #          wave,
 #          age_orig = age,
 #          age_test = age_imputed,
 #          state_test = state,
 #          int_date,
 #          birth_date,
 #          raddate,
 #          age_death,
 #          iwstat,
 #          iwmid) |>   
 #   right_join(hrs_jordan,
 #              by = join_by(hhidpn, wave)) |>    
 #   mutate(age_check = age_test - age_jordan,
 #          state_check = state_test - state_jordan) |> 
 #   filter(abs(zapsmall(state_check))>.2) 
   
# 
# 
# cog |> 
#    filter(hhidpn == "010099010") |> 
#   pivot_longer(-hhidpn,names_to = "year", values_to = "state") |> 
#   mutate(year = parse_number(year)) |> 
#   arrange(year)
 # early_dementia <- cog |>
 #   pivot_longer(cols = starts_with("cogfunction"),
 #                names_to = "year",
 #                names_pattern = "cogfunction(\\d{4})",
 #                values_to = "cogfunction") |>
 #   mutate(year = as.integer(year)) |>
 #   filter(cogfunction == 3) |>
 #   group_by(hhidpn) |>
 #   summarise(first_dementia_year = min(year), .groups = "drop") |>
 #   filter(first_dementia_year < 2000)
 # early_dementia
# ------------------------------------------------------------------- #
# first pass processing
# hrs_msm <- hrs_joined |>
hrs_msm <- hrs_joined |> 
  # pick age range to fit to, based on plot of support.
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
  mutate(obstype = ifelse(state_msm == 3, 3, 1)) |>
  # ------------------------------------------------------------------- #
  # RT: After cleaning we have to remove a solitary obs. one more time
  # 190 persons, they provide no information
  group_by(hhidpn) |> 
  # supposed to remove solitary observations.
  filter(n() > 1) |>
  ungroup()
# ------------------------------------------------------------------- #
# here we will work with NA values in the data
nas <- hrs_msm |>
  is.na() |>
  colSums()

# 360 misses in race are legit. 
# 54 in education too
# there are 25287 entries with unknown int_date
# these cases are also missing the health covariate info
# but they have a state info
# What should we do with such cases?
# TR: states were interpolated for these cases. To the extent that health covariates
# are to be treated as "ever" or "destined to / ever" variables, then we should 
# impute.
nas[nas > 0]
# ------------------------------------------------------------------- #
# Since all health conditions are "ever had" 
# We can obtain extra 28181 from 99887 missing cases
# by grouping and filling
hrs_msm <- hrs_msm |>
  group_by(hhidpn) |>
  fill(c(hypertension, diabetes, 
         heart_disease, stroke),  .direction = "downup") |>
  ungroup()
# ------------------------------------------------------------------- #
# I suggest using newly created obs_date instead of int_date
# hrs_msm |> 
#   filter(!is.na(int_date)) |>  
#   is.na() |> 
#   colSums()

# create obs_date code by TR
hrs_to_fit <- hrs_msm |>
  # a test condition to check pandemic effect on mort trend
  mutate(
    birth_date         = as_date(birth_date, origin = "1960-01-01"),
    birth_date_decimal = decimal_date(birth_date),
    obs_date           = birth_date_decimal + age
  ) |>
  filter(obs_date < (as_date("2019-dec-31") |> decimal_date())) |>
  # factor variables
  mutate(across(c(female, race, hispanic, education,
                  hypertension, diabetes, heart_disease,
                  stroke, ever_dementia), ~ as.factor(.)))
# ------------------------------------------------------------------- #
# end prep