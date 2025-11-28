
source("R/00_dependencies.R")
source("R/01_functions_extra.R")
# ------------- #
# model 1       #
# ------------- #

probs <- read_csv("Data/model1/probs.csv.gz")
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
         group_cols = c("replicate", "female", "period5"), 
         trans_col = "transition") |> 
  ungroup() |> 
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, period5, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model1/e50.csv.gz")

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
         group_cols = c("replicate", "female", "period5"), 
         trans_col = "transition") |> 
  ungroup() |> 
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, period5, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model1/e65.csv.gz")
# expectancies, age 65
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
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, period5, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model1/e85.csv.gz")

rm(probs);gc()

# ------------- #
# model 2       #
# ------------- #

probs <- read_csv("Data/model2/probs.csv.gz")
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
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, year, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model2/e50.csv.gz")

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
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, year, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model2/e65.csv.gz")

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
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, year, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model2/e85.csv.gz")
rm(probs);gc()


# ------------- #
# model 3       #
# ------------- #
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
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, year, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model3/e50.csv.gz")
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
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, year, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model3/e65.csv.gz")
gc()
# expectancies, age 65
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
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, year, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") |> 
  write_csv(file = "Data/model3/e85.csv.gz")
rm(probs);gc()

# end