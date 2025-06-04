
max_wave <- 16
hrs_file <- if_else(max_wave == 16, "randhrs1992_2022v1.dta","randhrs1992_2020v2.dta")
vars <- c("hhidpn", "hhid", "pn", "hacohort", "rabdate", "raddate",
          paste0("r", 5:max_wave, "iwstat"),
          paste0("r", 5:max_wave, "iwmid"))

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

hrs_in |> 
  pivot_longer(
  cols = matches("^r\\d{1,2}(iwstat|iwmid)$"),
  names_to = c("wave", ".value"),
  names_pattern = "r(\\d{1,2})([a-z]+)") |>
  mutate(wave = as.integer(wave)) |>
  filter(iwstat %in% c(1, 4, 5)) |> 
  mutate(iwmid = if_else(iwstat == 5, raddate, iwmid)) |> 
  group_by(hhidpn) |> 
  filter(n() > 1) |> 
  ungroup() |> 
  arrange(hhidpn,wave) |> 
  group_by(hhidpn) |> 
  mutate(wave_interval = lead(iwmid) - iwmid) |> 
  filter(!is.na(wave_interval)) |> 
  ggplot(aes(x = wave_interval)) +
  geom_density() +
  theme_minimal() +
  facet_wrap(~wave) +
  geom_vline(xintercept = 365.25*2,color="red")
