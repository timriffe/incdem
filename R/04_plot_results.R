# previous step
# source("R/01_data_preparation.R")
# source("R/02_estimate_models.R")
# source("R/03_predict_from_models.R")
# OR load prepared data
result_df <- read_csv("Data/results.csv.gz")
# ------------------------------------------------------------------- #
# plot rates
result_df |>
  filter(parse_number(to) > parse_number(from)) |> 
  mutate(Period = as.factor(int_date_decimal),
         transition = paste0("haz",parse_number(from),parse_number(to))) |>
  ggplot(aes(x        = age, 
             y        = rate, 
             color    = transition, 
             linetype = Period)) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    subtitle = "*Linear time trend",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap(~sex) + 
  theme(legend.position = "bottom")
# ------------------------------------------------------------------- #
# plot probability
result_df |>
  mutate(Period = as.factor(int_date_decimal)) |> 
  ggplot(aes(x        = age, 
             y        = prob, 
             color    = to, 
             linetype = Period)) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    subtitle = "*Linear time trend",
    x        = "Age",
    y        = "Hazard Rate",
    color   = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() +
  facet_wrap(sex ~ from) + 
  theme(legend.position = "bottom")
# ------------------------------------------------------------------- #
# hazard rate
result_df |>
  mutate(
    from = from |> as.character() |> parse_number(),
    to   =  to  |> as.character() |> parse_number(),
    transition = paste(from, "→", to)
  ) |>
  filter(to   == 2,
         from == 1) |>
  ggplot(aes(
    x        = age,
    y        = rate,
    color    = as.factor(int_date_decimal),
    linetype = sex
  )) +
  geom_line(linewidth = 1) +
  labs(
    title    = "Estimated Transition Hazards by Age",
    subtitle = "*Linear time trend",
    x        = "Age",
    y        = "Hazard Rate",
    color    = "Transition"
  ) +
  theme_minimal() +
  scale_y_log10() + 
  theme(legend.position = "bottom")
# ------------------------------------------------------------------- #
# check support; this was used to determine which age range
# to fit over. Reason: spline tail misbehaves in young ages
# if all ages are considered.
hrs_msm |>
  filter(state_clean == 2) |>
  count(age_int = floor(age)) |>
  ggplot(aes(x = age_int, 
             y = n)) +
  geom_col() +
  labs(title = "Number of Observations in Dementia by Age",
       x     = "Age", 
       y     = "Count")
# ------------------------------------------------------------------- #