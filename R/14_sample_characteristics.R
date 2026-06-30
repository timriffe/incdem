# Load required packages
source("R/00_dependencies.R")
source("R/01_functions.R")
# Select necessary variables from main dataset
max_wave <- 16
hrs_file <- if_else(max_wave == 16, 
                    "randhrs1992_2022v1.dta",
                    "randhrs1992_2020v2.dta")
hrs_vars <- c("hhidpn", "hhid", "pn", 
              "hacohort", "rabdate", "raddate",
              paste0("r", 1:max_wave, "iwstat"),
              paste0("r", 1:max_wave, "iwmid"),
              paste0("r", 1:max_wave, "hibpe"),
              paste0("r", 1:max_wave, "diabe"),
              paste0("r", 1:max_wave, "hearte"),
              paste0("r", 1:max_wave, "stroke"),
              paste0("r", 1:max_wave, "wtcrnh"),
              "ragender", "raracem", "rahispan", "rabyear", "raeduc")

# Load main RAND HRS dataset
hrs_in <- read_dta(file.path("Data",hrs_file), 
                   encoding   = "UTF-8", 
                   col_select = all_of(hrs_vars)) |> 
  zap_labels() |> 
  mutate(
    hhid   = as.character(hhid),
    pn     = as.character(pn),
    hhidpn = paste0(hhid, pn)
  )


# Load cognition dataset
cog_colnames <- names(read_dta("Data/cogfinalimp_9520wide.dta", n_max = 0))
cog_cols_needed <- c("hhid", "pn", grep("^cogfunction\\d{4}$", 
                                        cog_colnames, 
                                        value = TRUE))

cog <- read_dta("./Data/cogfinalimp_9520wide.dta", 
                encoding   = "UTF-8", 
                col_select = all_of(cog_cols_needed)) |>
  zap_labels() |>
  mutate(
    hhid   = as.character(hhid),
    pn     = as.character(pn),
    hhidpn = paste0(hhid, pn)
  ) |>
  select(hhidpn, starts_with("cogfunction"))

# Pre-process HRS dataset
hrs_long <- hrs_in |>
  mutate(
    female = as.integer(ragender == 2),
    race   = case_when(
      raracem == 1 ~ "white",
      raracem == 2 ~ "black",
      raracem == 3 ~ "other"
    ),
    hispanic = if_else(rahispan == 1, "hispanic", "non-hispanic"),
    education = case_when(
      raeduc %in% c(1, 2) ~ "ls hs/ged",
      raeduc == 3         ~ "hs",
      raeduc %in% c(4, 5) ~ ">hs"
    ),
    birth_year = rabyear,
    cohort     = hacohort
  ) |>
  pivot_longer(
    cols = matches("^r\\d{1,2}(iwstat|iwmid|hibpe|diabe|hearte|stroke|wtcrnh)$"),
    names_to      = c("wave", ".value"),
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
  mutate(n = n()) |> 
  ungroup() |> 
  filter(n > 1) |> 
  arrange(hhidpn, wave)

z <- hrs_long %>%
  group_by(hhidpn) |>
  mutate(age_imputed = impute_age(age, wave)) %>% 
  select(hhidpn, age, age_imputed)

# % of ages missing in original sample that were imputed by impute_age
(z$age %>% is.na() %>% sum()) / length(z$age)


# ------------------------------------------------------------------- #
# Here I calculate sample characteristics
rm(list = ls())
source("R/00_dependencies.R")
source("R/01_functions.R")
source("R/02_prepare_hrs.R")
# sample characteristics
# number of unique persons
hrs_to_fit_prepped <- hrs_to_fit |>
  mutate(period5 = case_when(between(obs_date, 2004, 2010) ~ "period 1",
                             between(obs_date, 2010, 2015) ~ "period 2",
                             obs_date > 2015 ~ "period 3",
                             TRUE ~ "other"),
         year = floor(obs_date)) |>
  select(wave, 
         hhidpn, female, education,
         pwt = wtcrnh, age, hypertension:stroke,
         state_msm, ever_dementia, obstype, period5,
         year) |> 
  filter(period5 != "other",
         year > 2003,
         state_msm != 3) |> 
  mutate(period = as.factor(period5))

# number of people in sample 
hrs_to_fit_prepped %>% 
  pull(hhidpn) %>% 
  unique() %>%
  length()

# sex distribution
hrs_to_fit_prepped %>% 
  select(hhidpn, female) %>%
  distinct() %>% 
  janitor::tabyl(female)

# mean sd age at baseline (first appearance)
hrs_to_fit_prepped %>% 
  select(hhidpn, female, age) %>%
  filter(age == min(age), .by = "hhidpn") %>%
  summarise(m_age = mean(age),
            s_age = sd(age))

# by sex
hrs_to_fit_prepped %>% 
  select(hhidpn, female, age) %>%
  filter(age == min(age), .by = "hhidpn") %>%
  summarise(m_age = mean(age),
            s_age = sd(age), .by = "female")

# distribution of individual follow-up duration (years) Exposure?
exposure <- hrs_to_fit_prepped %>% 
  select(hhidpn, age) %>%
  group_by(hhidpn) %>%
  summarise(
    df = max(age) - min(age),
    .groups = "drop"
  ) %>% 
  filter(df > 0)

mean(exposure$df)
sd(exposure$df)

# wave spacing
wave_spacing <- hrs_to_fit_prepped %>% 
  select(hhidpn, year) %>%
  arrange(hhidpn, year) %>%
  group_by(hhidpn) %>%
  mutate(
    spacing = year - lag(year)
  ) %>%
  ungroup() %>%
  filter(!is.na(spacing)) %>% 
  filter(spacing > 0)

janitor::tabyl(wave_spacing$spacing)
