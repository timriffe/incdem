

# using a given model, let's calculate the hazard ratios of
# 1. change 2019/2004
# 2. sex (m/f)

# To do this, we need to age-standardize somehow. It makes sense to me (TR) to age standardize using the corresponding stationary occupancy times. That is, use lh to weight HD and HU, and use lu to weight UD. Makes sense to me anyway. For time, we can take the avg of ls and the start and end. For sex, we can take the avg as well, and we'll do it for just 2004 and 2019 to keep it simple.

source("R/00_dependencies.R")
source("R/01_functions_extra.R")

model <- 2
from_age <- 50
age_interval <- 0.25
years_keep <- c(2004, 2019)

data_dir <- paste0("Data/model", model)

# transitions:
# 1 -> 2: HU, standardized with Lx1 / Lh
# 1 -> 3: HD, standardized with Lx1 / Lh
# 2 -> 3: UD, standardized with Lx2 / Lu

transitions_keep <- tibble::tribble(
  ~from, ~to, ~transition, ~standard_state,
  1,     2,   "HU",        "Lx1",
  1,     3,   "HD",        "Lx1",
  2,     3,   "UD",        "Lx2"
)

# ------------- #
# read and trim #
# ------------- #

haz <- readr::read_csv(file.path(data_dir, "adj_haz_replicates.csv.gz"), show_col_types = FALSE) |>
  dplyr::filter(
    year %in% years_keep,
    age >= from_age
  ) |>
  dplyr::inner_join(
    transitions_keep,
    by = c("from", "to")
  )
gc()

probs <- readr::read_csv(file.path(data_dir, "probs.csv.gz"), show_col_types = FALSE) |>
  dplyr::filter(
    year %in% years_keep,
    age >= from_age
  )
gc()
# ------------- #
# Lx standards  #
# ------------- #

lxs <- calc_lxs(
  probs = probs,
  from_age = from_age,
  age_interval = age_interval,
  init = c(`1` = 1, `2` = 0),
  init_method = "init",
  from_col = "from",
  to_col = "to",
  age_col = "age",
  p_col = "p",
  group_cols = c("replicate", "female", "year"),
  trans_col = "transition"
)

rm(probs)
gc()

lxs_long <- lxs |>
  dplyr::select(replicate, female, year, age, Lx1, Lx2) |>
  tidyr::pivot_longer(
    cols = c(Lx1, Lx2),
    names_to = "standard_state",
    values_to = "Lx"
  )

rm(lxs)
gc()

# ------------- #
# change HRs    #
# 2019 / 2004   #
# age standard: #
# same within   #
# sex and rep   #
# ------------- #

change_standards <- lxs_long |>
  dplyr::group_by(replicate, female, age, standard_state) |>
  dplyr::summarise(
    Lx_std = mean(Lx),
    .groups = "drop"
  )

haz_std <- haz |>
  dplyr::left_join(
    change_standards,
    by = c("replicate", "female", "age", "standard_state")
  ) |>
  dplyr::group_by(replicate, female, year, transition) |>
  dplyr::summarise(
    haz_std = sum(haz * Lx_std, na.rm = TRUE) / sum(Lx_std, na.rm = TRUE),
    .groups = "drop"
  )

change_hrs <- haz_std |>
  tidyr::pivot_wider(
    names_from = year,
    values_from = haz_std,
    names_prefix = "haz_"
  ) |>
  dplyr::mutate(
    HR_2019_2004 = haz_2019 / haz_2004
  ) |>
  dplyr::select(
    replicate, female, transition,
    haz_2004, haz_2019, HR_2019_2004
  )

change_summary <- change_hrs |>
  dplyr::group_by(female, transition) |>
  dplyr::summarise(
    HR = stats::median(HR_2019_2004, na.rm = TRUE),
    lower = stats::quantile(HR_2019_2004, 0.025, na.rm = TRUE),
    upper = stats::quantile(HR_2019_2004, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# ------------- #
# sex HRs       #
# women / men   #
# age standard: #
# same within   #
# year and rep  #
# ------------- #

sex_standards <- lxs_long |>
  dplyr::group_by(replicate, year, age, standard_state) |>
  dplyr::summarise(
    Lx_std = mean(Lx),
    .groups = "drop"
  )

sex_rates <- haz |>
  dplyr::left_join(
    sex_standards,
    by = c("replicate", "year", "age", "standard_state")
  ) |>
  dplyr::group_by(replicate, year, female, transition) |>
  dplyr::summarise(
    haz_std = sum(haz * Lx_std, na.rm = TRUE) / sum(Lx_std, na.rm = TRUE),
    .groups = "drop"
  )

sex_hrs <- sex_rates |>
  tidyr::pivot_wider(
    names_from = female,
    values_from = haz_std,
    names_prefix = "haz_female_"
  ) |>
  dplyr::mutate(
    HR_female_male = haz_female_0 / haz_female_1
  ) |>
  dplyr::select(
    replicate, year, transition,
    haz_female_0, haz_female_1, HR_female_male
  )

sex_summary <- sex_hrs |>
  dplyr::group_by(year, transition) |>
  dplyr::summarise(
    HR = stats::median(HR_female_male, na.rm = TRUE),
    lower = stats::quantile(HR_female_male, 0.025, na.rm = TRUE),
    upper = stats::quantile(HR_female_male, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


rm(
  haz,
  lxs_long,
  change_standards,
  change_rates,
  change_hrs,
  change_summary,
  sex_standards,
  sex_rates,
  sex_hrs,
  sex_summary
)

gc()

# end
