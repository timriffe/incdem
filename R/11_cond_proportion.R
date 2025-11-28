source("R/00_dependencies.R")
source("R/01_functions.R")
source("R/01_functions_extra.R")
# Proportion of life lived with dementia at age 65 and 85 (a simple, policy-friendly metric)
# TR: I interpret this as "upward-looking" at ages 65 and 85:

# -------------------- #
# model 1              #
# -------------------- #

probs <- read_csv("Data/model1/probs.csv.gz")
# proportion healthy, age 50
calc_exs(probs = probs,
         from_age = 50, 
         age_interval = 0.25, 
         init = c(`1` = 1, `2` = 0) ,
         init_method = "init", 
         from_col = "from", 
         to_col = "to", 
         age_col = "age", 
         p_col = "p", 
         group_cols = c("replicate", "female", "period5"), 
         trans_col = "transition") |> 
  ungroup() |> 
  mutate(prop = DFLE / (DFLE + DLE)) |> 
  select(-DFLE,-DLE) |> 
  group_by(female, period5) |> 
  summarize(prop50 = median(prop),
            lower = quantile(prop,0.025),
            upper = quantile(prop,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model1/prop50.csv.gz")
gc()
#  proportion healthy, age 65
calc_exs(probs = probs,
         from_age = 65, 
         age_interval = 0.25, 
         init = c(`1` = 1, `2` = 0) ,
         init_method = "init", 
         from_col = "from", 
         to_col = "to", 
         age_col = "age", 
         p_col = "p", 
         group_cols = c("replicate", "female", "period5"), 
         trans_col = "transition") |> 
  ungroup() |> 
  mutate(prop = DFLE / (DFLE + DLE)) |> 
  select(-DFLE,-DLE) |> 
  group_by(female, period5) |> 
  summarize(prop50 = median(prop),
            lower = quantile(prop,0.025),
            upper = quantile(prop,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model1/prop65.csv.gz")
gc()
#  proportion healthy, age 85
calc_exs(probs = probs,
         from_age = 85, 
         age_interval = 0.25, 
         init = c(`1` = 1, `2` = 0) ,
         init_method = "init", 
         from_col = "from", 
         to_col = "to", 
         age_col = "age", 
         p_col = "p", 
         group_cols = c("replicate", "female", "period5"), 
         trans_col = "transition") |> 
  ungroup() |> 
  mutate(prop = DFLE / (DFLE + DLE)) |> 
  select(-DFLE,-DLE) |> 
  group_by(female, period5) |> 
  summarize(prop50 = median(prop),
            lower = quantile(prop,0.025),
            upper = quantile(prop,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model1/prop85.csv.gz")
gc()
rm(probs);gc()
# -------------------- #
# model 2              #
# -------------------- #
probs <- read_csv("Data/model2/probs.csv.gz")
#  proportion healthy, age 50
calc_exs(probs = probs,
         from_age = 50, 
         age_interval = 0.25, 
         init = c(`1` = 1, `2` = 0) ,
         init_method = "init", 
         from_col = "from", 
         to_col = "to", 
         age_col = "age", 
         p_col = "p", 
         group_cols = c("replicate", "female", "year"), 
         trans_col = "transition") |> 
  ungroup() |> 
  mutate(prop = DFLE / (DFLE + DLE)) |> 
  select(-DFLE,-DLE) |> 
  group_by(female, year) |> 
  summarize(prop50 = median(prop),
            lower = quantile(prop,0.025),
            upper = quantile(prop,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model2/prop50.csv.gz")
gc()
# expectancies, age 65
calc_exs(probs = probs,
         from_age = 65, 
         age_interval = 0.25, 
         init = c(`1` = 1, `2` = 0) ,
         init_method = "init", 
         from_col = "from", 
         to_col = "to", 
         age_col = "age", 
         p_col = "p", 
         group_cols = c("replicate", "female", "year"), 
         trans_col = "transition") |> 
  ungroup() |> 
  mutate(prop = DFLE / (DFLE + DLE)) |> 
  select(-DFLE,-DLE) |> 
  group_by(female, year) |> 
  summarize(prop50 = median(prop),
            lower = quantile(prop,0.025),
            upper = quantile(prop,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model2/prop65.csv.gz")
gc()
# expectancies, age 85
calc_exs(probs = probs,
         from_age = 85, 
         age_interval = 0.25, 
         init = c(`1` = 1, `2` = 0) ,
         init_method = "init", 
         from_col = "from", 
         to_col = "to", 
         age_col = "age", 
         p_col = "p", 
         group_cols = c("replicate", "female", "year"), 
         trans_col = "transition") |> 
  ungroup() |> 
  mutate(prop = DFLE / (DFLE + DLE)) |> 
  select(-DFLE,-DLE) |> 
  group_by(female, year) |> 
  summarize(prop50 = median(prop),
            lower = quantile(prop,0.025),
            upper = quantile(prop,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model2/prop85.csv.gz")

rm(probs);gc()
# -------------------- #
# model 3              #
# -------------------- #

probs <- read_csv("Data/model3/probs.csv.gz")
# expectancies, age 50
calc_exs(probs = probs,
         from_age = 50, 
         age_interval = 0.25, 
         init = c(`1` = 1, `2` = 0) ,
         init_method = "init", 
         from_col = "from", 
         to_col = "to", 
         age_col = "age", 
         p_col = "p", 
         group_cols = c("replicate", "female", "year"), 
         trans_col = "transition") |> 
  ungroup() |> 
  mutate(prop = DFLE / (DFLE + DLE)) |> 
  select(-DFLE,-DLE) |> 
  group_by(female, year) |> 
  summarize(prop50 = median(prop),
            lower = quantile(prop,0.025),
            upper = quantile(prop,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model3/prop50.csv.gz")
gc()
# expectancies, age 65
calc_exs(probs = probs,
         from_age = 65, 
         age_interval = 0.25, 
         init = c(`1` = 1, `2` = 0) ,
         init_method = "init", 
         from_col = "from", 
         to_col = "to", 
         age_col = "age", 
         p_col = "p", 
         group_cols = c("replicate", "female", "year"), 
         trans_col = "transition") |> 
  ungroup() |> 
  mutate(prop = DFLE / (DFLE + DLE)) |> 
  select(-DFLE,-DLE) |> 
  group_by(female, year) |> 
  summarize(prop50 = median(prop),
            lower = quantile(prop,0.025),
            upper = quantile(prop,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model3/prop65.csv.gz")
gc()
# expectancies, age 85
calc_exs(probs = probs,
         from_age = 85, 
         age_interval = 0.25, 
         init = c(`1` = 1, `2` = 0) ,
         init_method = "init", 
         from_col = "from", 
         to_col = "to", 
         age_col = "age", 
         p_col = "p", 
         group_cols = c("replicate", "female", "year"), 
         trans_col = "transition") |> 
  ungroup() |> 
  mutate(prop = DFLE / (DFLE + DLE)) |> 
  select(-DFLE,-DLE) |> 
  group_by(female, year) |> 
  summarize(prop50 = median(prop),
            lower = quantile(prop,0.025),
            upper = quantile(prop,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model3/prop85.csv.gz")

rm(probs);gc()

# end