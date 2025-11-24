# previous step
source("R/00_package_and_functions.R")

# ----------------------- #
# model 1                 #
# ----------------------- #
probs <- read_csv("Data/model1/probs.csv.gz")
# ------------------------------------------------------------------- #

probs |> 
  filter(to > from) |> 
  group_by(period5, female, age, from, to) |> 
  summarize(p_med = median(p),
            lower = quantile(p, 0.025),
            upper = quantile(p, 0.975)) |> 
  ggplot(aes(x = age, y = p_med, color = interaction(from, to))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = interaction(from, to)), color = "transparent", alpha = .3) +
  facet_wrap(female~period5) +
  scale_y_log10()

# ----------------------- #
# model 2                 #
# ----------------------- #
probs <- read_csv("Data/model2/probs.csv.gz")
# ------------------------------------------------------------------- #
probs |> 
  filter(to > from) |> 
  group_by(female, year, age, from, to) |> 
  summarize(p_med = median(p),
            lower = quantile(p, 0.025),
            upper = quantile(p, 0.975)) |> 
  ggplot(aes(x = age, y = p_med, color = interaction(from, to))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = interaction(from, to)), color = "transparent", alpha = .3) +
  facet_wrap(year~female) +
  scale_y_log10()

rm(probs);gc()

# ----------------------- #
# model 3                 #
# ----------------------- #
probs <- read_csv("Data/model3/probs.csv.gz")
# ------------------------------------------------------------------- #
probs |> 
  filter(to > from) |> 
  group_by(female, year, age, from, to) |> 
  summarize(p_med = median(p),
            lower = quantile(p, 0.025),
            upper = quantile(p, 0.975)) |> 
  ggplot(aes(x = age, y = p_med, color = interaction(from, to))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = interaction(from, to)), color = "transparent", alpha = .3) +
  facet_wrap(year~female) +
  scale_y_log10()

rm(probs);gc()