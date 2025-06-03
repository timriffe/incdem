# Load required packages
source("R/00_packages.R")

impute_age <- function(age, wave){
  
  age <-  zoo::na.approx(age, x = wave, na.rm = FALSE)
  if (any(is.na(age))){
    if (is.na(age[1])){
      # handle leading NAs
      diffs            <- diff(wave) * 2
      nNAs             <- rle(is.na(age))$lengths[1]
      first_non_NA_age <- age[!is.na(age)][1]
      subtract_this    <- rev(1:nNAs * diffs[1:nNAs]) 
      age[1:nNAs]      <- first_non_NA_age - subtract_this
    }
    n <- length(age)
    if (is.na(age[n])){
      diffs            <- diff(wave) * 2
      nNAs             <- rle(is.na(rev(age)))$lengths[1]
      last_non_NA_age  <- age[!is.na(age)] |> rev() %>% '['(1)
      add_this         <- (rev(diffs)[1:nNAs] |> rev()) * 1:nNAs
      age[(n-nNAs+1):n]  <- last_non_NA_age + add_this
    }
  }
  age
}
# Select necessary variables from main dataset
vars <- c("hhidpn", "hhid", "pn", "hacohort", "rabdate", "raddate",
          paste0("r", 5:16, "iwstat"),
          paste0("r", 5:16, "iwmid"),
          paste0("r", 5:16, "hibpe"),
          paste0("r", 5:16, "diabe"),
          paste0("r", 5:16, "hearte"),
          paste0("r", 5:16, "stroke"),
          "ragender", "raracem", "rahispan", "rabyear", "raeduc")
# hrs_file <- "randhrs1992_2020v2.dta"
hrs_file <- "randhrs1992_2022v1.dta"
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
  arrange(hhidpn, wave) |>
  group_by(hhidpn) |>
  mutate(age_imputed = impute_age(age, wave)) |> 
  ungroup() |>
   group_by(hhidpn) |>
  mutate(across(c(hibpe, diabe, hearte, stroke),
                ~ cummax(replace_na(., 0)),
                .names = "{.col}_status")) |>
  ungroup()


# stata code:
# * recode 1/2 = 0 (no dementia), 3 = 1 (dementia)
# recode cogfunction (1/2 = 0) (3 = 1), gen(cog_dementia)
# 
# * forward fill missing values
# bysort hhid pn (wave): replace cog_dementia = cog_dementia[_n-1] if missing(cog_dementia)
# 
# * backward fill remaining missing values
# gsort hhid pn -wave
# by hhid pn: replace cog_dementia = cog_dementia[_n-1] if missing(cog_dementia)
# sort hhid pn wave
# 
# * generate running maximum to preserve dementia status
# gen byte dementia = cog_dementia
# by hhid pn (wave): replace dementia = cond(dementia[_n-1] > cog_dementia, dementia[_n-1], cog_dementia) if _n > 1

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
  rename(birth_date = rabdate)

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
 #          age_death,
 #          iwstat) |>   
 #   right_join(hrs_jordan,
 #              by = join_by(hhidpn, wave)) |>    
 #   mutate(age_check = age_test - age_jordan) |>    
 #   # filter(abs(zapsmall(age_check))>.2) 
 #   filter(!is.na(age_orig), abs(age_check)>.1) 
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
