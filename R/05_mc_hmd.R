source("R/00_package_and_functions.R")
library(HMDHFDplus)
library(janitor)

Deaths <- readHMDweb("USA","Deaths_1x1", 
                 username = Sys.getenv("us"), 
                 password = Sys.getenv("pw")) |> 
  clean_names() |> 
  select(-total, -open_interval) |> 
  filter(between(age, 40, 110),   # these could be global variables set elsewhere
         between(year, 2004, 2019)) |> 
  pivot_longer(c(male, female),names_to = "sex", values_to = "deaths")

Exposures <- readHMDweb("USA","Exposures_1x1", 
                 username = Sys.getenv("us"), 
                 password = Sys.getenv("pw"))|> 
  clean_names() |> 
  select(-total, -open_interval) |> 
  filter(between(age, 40, 110),   # these could be global variables set elsewhere
         between(year, 2004, 2019)) |> 
  pivot_longer(c(male, female), names_to = "sex", values_to = "exposure")


# make a helper that smooths mortality ages 40-110,
# interpolating to .25 age intervals (.25 should also be a global parameter)

# when done, select down to ages 50-100



external_ref <-
  Deaths |> 
  left_join(Exposures, by = join_by(year, age, sex)) |> 
  mutate(female = as.factor(as.integer(sex == "female"))) |> 
  select(-sex) |> 
  filter(between(age, 40, 110)) |>
  mutate(deaths = round(deaths)) 

gam_chunk <- function(df, m = 1, k = 11, age_interval_out = 0.25) {

  age_min <- min(df$age, na.rm = TRUE)
  age_max <- max(df$age, na.rm = TRUE)
  
  # data for fitting
  data_inside <- df %>%
    select(age, deaths, exposure, year)
  
  # prediction grid
  age_grid <- expand_grid(
    age      = seq(age_min, age_max, by = age_interval_out),
    year     = sort(unique(df$year)),
    exposure = 1
  )
  
  fit <- mgcv::gam(
    deaths ~ s(age, bs = "ps", m = m, k = k) + year,
    offset = log(exposure),
    family = poisson(link = "log"),
    data   = data_inside
  )
  
  mu_hat <- predict(fit, newdata = age_grid, type = "response")
  
  # df-out: one row per (age, year) with smoothed deaths (or rates)
  age_grid %>%
    mutate(mx = mu_hat)
}

external_ref_quarter <- external_ref %>%
  group_by(female) %>%
  group_modify(~ gam_chunk(.x, 
                           m = 1, 
                           k = 11, 
                           age_interval_out = 0.25)) %>%
  ungroup()

external_ref |> 
  mutate(mx = deaths / exposure) |> 
ggplot(aes(x=age,y=mx,color=female)) +
  geom_point(size = .4) +
  geom_line(data = external_ref_quarter) +
  facet_wrap(~year) +
  scale_y_log10()

nrow(external_ref)
N <- 1000

# Create Monte Carlo ensemble
external_ref_mc <-
  tibble(replicate = 1:N) |> 
  cross_join(external_ref) |> 
  arrange(replicate,female, year, age) |> 
  rename(deaths_ref = deaths) |> 
  mutate(deaths = rpois(n(),deaths_ref)) |> 
  group_by(replicate, female) |> 
  group_modify(~ gam_chunk(.x, 
                           m = 1, 
                           k = 11, 
                           age_interval_out = 0.25)) %>%
  ungroup() |> 
  filter(between(age,50,100)) |> 
  mutate(mx = as.vector(mx))
write_csv(external_ref_mc, file = "Data/external_ref_mc.csv.gz")

rm(external_ref_mc);rm(external_ref_quarter);gc()
# VERY tight bounds
# external_ref_mc |> 
#   filter(year == 2015) |> 
#   group_by(year,female,age) |> 
#   summarize(mx_med = median(mx),
#             lower = quantile(mx,0.025),
#             upper = quantile(mx,0.975)) |> head()
#   ggplot(aes(x = age,y=mx_med, color = female)) +
#   geom_line() +
#   scale_y_log10() +
#   geom_ribbon(aes(ymin=lower,ymax=upper,fill = female),alpha=.3,color="transparent") 


