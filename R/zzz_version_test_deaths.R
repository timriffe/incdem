
source("R/00_packages.R")


vars <- c("hhid", "pn","rabdate","raddate")

hrs_1 <- read_dta(file.path("Data","randhrs1992_2020v2.dta"), 
                  encoding = "UTF-8", 
                  col_select = all_of(vars)) |> 
  zap_labels() |> 
  mutate(
    hhid = as.character(hhid),
    pn = as.character(pn),
    hhidpn = paste0(hhid, pn)
  )

hrs_2 <- read_dta(file.path("Data","randhrs1992_2022v1.dta"), 
                   encoding = "UTF-8", 
                   col_select = all_of(vars)) |> 
  zap_labels() |> 
  mutate(
    hhid = as.character(hhid),
    pn = as.character(pn),
    hhidpn = paste0(hhid, pn)
  )

deaths1 <-
  hrs_1 |> 
  filter(!is.na(raddate)) |> 
  mutate(death_date = as_date(raddate, origin = "1960-01-01"),
         birth_date = as_date(rabdate, origin = "1960-01-01"),
         death_year = year(death_date)) |> 
  group_by(death_year) |> 
  summarize(deaths = n()) |> 
  mutate(version = "randhrs1992_2020v2")

deaths2 <-
  hrs_2 |> 
  filter(!is.na(raddate)) |> 
  mutate(death_date = as_date(raddate, origin = "1960-01-01"),
         birth_date = as_date(rabdate, origin = "1960-01-01"),
         death_year = year(death_date)) |> 
  group_by(death_year) |> 
  summarize(deaths = n()) |> 
  mutate(version = "randhrs1992_2022v1")

bind_rows(deaths1,deaths2) |> 
  ggplot(aes(x = death_year, y = deaths, color = version)) +
  geom_line() +
  geom_vline(xintercept = 2019)
