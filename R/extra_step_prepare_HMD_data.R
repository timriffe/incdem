# -----------------------------------------------------------------------------#
source("R/00_package_and_functions.R")

# -----------------------------------------------------------------------------#
# load HMD data Death counts
# nothing interesting happens here just data wrangling
Dx <- read_table("Data/Deaths_1x1.txt", skip = 1) |>
  dplyr::select(-Total) |>
  mutate(age = parse_number(Age)) |>
  filter(Year %in% c(2000:2020))  |>
  mutate(period = case_when(between(Year, 2000, 2006) ~ "period 1",
                            between(Year, 2006, 2012) ~ "period 2",
                            between(Year, 2012, 2020) ~ "period 3"),
         period = as.factor(period)) |>
  group_by(period, age) |>
  summarise(across(c(Female, Male), ~ sum(.)), .groups = "drop") |>
  dplyr::select(period, age, Male, Female) |>
  pivot_longer(-c(period, age),
               names_to  = "female",
               values_to = "Dx") |>
  mutate(female = ifelse(female == "Female", 1, 0))

# -----------------------------------------------------------------------------#
# HMD Exposures
Ex <- read_table("Data/Exposures_1x1.txt", skip = 1) |>
  dplyr::select(-Total) |>
  mutate(age = parse_number(Age)) |>
  filter(Year %in% c(2000:2020))  |>
  mutate(period = case_when(between(Year, 2000, 2006) ~ "period 1",
                            between(Year, 2006, 2012) ~ "period 2",
                            between(Year, 2012, 2020) ~ "period 3"),
         period = as.factor(period)) |>
  group_by(period, age) |>
  summarise(across(c(Female, Male), ~ sum(.)), .groups = "drop") |>
  dplyr::select(period, age, Male, Female) |>
  pivot_longer(-c(period, age),
               names_to  = "female",
               values_to = "Ex") |>
  mutate(female = ifelse(female == "Female", 1, 0))

# -----------------------------------------------------------------------------#
# calculate mx
hmd_mx <- Dx |>
  full_join(Ex) |>
  mutate(mu  = Dx / Ex)

# -----------------------------------------------------------------------------#
# we do periods, so no need for HMD LT
# # hmd lifetable
# hmd_m <- read_table("Data/mltper_1x1.txt", skip = 1) |>
#   filter(Year %in% c(2000:2020)) |>
#   mutate(period = case_when(between(Year, 2000, 2006) ~ "period 1",
#                             between(Year, 2006, 2012) ~ "period 2",
#                             between(Year, 2012, 2020) ~ "period 3"),
#          period = as.factor(period)) |>
#   mutate(age = parse_number(Age)) |>
#   dplyr::select(period, age, mu_hmd = mx) |>
#   mutate(female = 0)
# 
# hmd_f <- read_table("Data/fltper_1x1.txt", skip = 1) |>
#   filter(Year %in% c(2000:2020)) |>
#   mutate(period = case_when(between(Year, 2000, 2006) ~ "period 1",
#                             between(Year, 2006, 2012) ~ "period 2",
#                             between(Year, 2012, 2020) ~ "period 3"),
#          period = as.factor(period)) |>
#   mutate(age = parse_number(Age)) |>
#   dplyr::select(period, age, mu_hmd = mx) |>
#   mutate(female = 1)

# full lt
# hmd_lt <- hmd_f %>% 
#   full_join(hmd_m)

# Merge 
# hmd <- hmd_mx %>% 
#   full_join(hmd_lt)
# -----------------------------------------------------------------------------#
# smoothing data
# NOTE THIS PART IS JUST A TEST
# 0.25 age smoothing comes later line 144
hmd <- hmd_mx |>
  # keep ages 50-100
  filter(between(age, 50, 100)) |>
  group_by(female) |>
  # since we are using poisson family pclm we need tp dound Deaths
  mutate(Dx = round(Dx)) |>
  group_nest() |>
  mutate(result = map(data, ~
      gam(
        Dx ~ s(age, bs = "ps", m = 1, k = 11) + period,
        offset = log(Ex),
        family = poisson(link = "log"),
        data = .x
      )),
  data = map2(data, result, ~ .x |>
                mutate(mu_hat = predict(
                  .y, newdata = .x, type = "response"
                )))) |>
  dplyr::select(-result) |>
  unnest(data)

# looks good to me
# NOTE: We can further improve fit by using using independent periods
# but I think passing it as a covariate is more motivated demographically
hmd %>%
  ggplot(aes(x = age)) +
  geom_line(aes(
    y = mu,
    color = "Observed",
    linetype = "Observed"
  ), linewidth = 1) +    # Observed
  # geom_line(aes(
  #   y = mu_hmd,
  #   color = "HMD",
  #   linetype = "HMD"
  # ), linewidth = 1) + # HMD
  geom_line(aes(
    y = mu_hat,
    color = "Smooth",
    linetype = "Smooth"
  ), linewidth = 1) +  # Smooth
  scale_y_log10() +
  labs(y = "Mortality Rate (Log Scale)", x = "Age") +
  scale_color_manual(values = c(
    "Observed"  = "black",
    # "HMD"       = "blue",
    "Smooth" = "red"
  )) +
  scale_linetype_manual(values = c(
    "Observed" = "solid",
    # "HMD"   = "dashed",
    "Smooth" = "dotted"
  )) +
  facet_wrap(period ~ female) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.key = element_rect(fill = "white", color = "white")
  ) +
  guides(color    = guide_legend(title = "Model"),
         linetype = guide_legend(title = "Model"))
# ----------------------------------------------------------------------- #
# graduate data to 0.25
# THIS IS THE RESULT
new_data <- expand_grid(age    = seq(50, 100, 0.25),
                        female = unique(hmd_mx$female),
                        period = unique(hmd_mx$period),
                        Ex     = 1)

# smooth and graduate data at the same time
# same model
hmd_25 <- hmd_mx |>
  ungroup() |> 
  filter(between(age, 50, 100)) |>
  mutate(Dx = round(Dx)) |>
  group_nest(female) |>
  nest_join(new_data, by = c("female")) |> 
  mutate(model = map(data, ~
                        gam(
                          Dx ~ s(age, bs = "ps", m = 1, k = 11) + period,
                          offset = log(Ex),
                          family = poisson(link = "log"),
                          data = .x
                        ))) |>
  mutate(pred = map2(.x = model, 
                     .y = new_data, ~ .y |>
                       mutate(new = predict(.x, newdata = .y, type = "response")))) %>% 
  dplyr::select(-c(data:model))  |>
  unnest(pred)  |>
  select(-Ex) |>
  rename(mu = new)  |>
  # calculate qx
  mutate(qx = 0.25 * mu  / (1 + 0.125 * mu))

# qx looks ok
hmd_25 |>
  geom_line(aes(x = age, y = qx, color = period)) +
  facet_wrap(~ female) +
  theme_bw() + 
  scale_y_log10()

# mx same
hmd_25 |>
  ggplot() +
  geom_line(aes(x = age, y = mu, color = period)) +
  facet_wrap(~ female) +
  theme_bw() + 
  scale_y_log10()
# save
save(hmd_25, file = "Data/hmd_mort.RData")

