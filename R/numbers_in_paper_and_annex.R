ker <- hrs_in |>
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


ker1 <- ker |>
  group_by(hhidpn) |>
  mutate(age_imputed = impute_age(age, wave))



ker1 %>% 
  filter(is.na(age_imputed))

ker1$age_imputed %>% is.na() %>% sum()

(ker$age %>% is.na() %>% sum()) / length(ker$age)

