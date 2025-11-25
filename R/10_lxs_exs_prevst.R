
source("R/00_dependencies.R")
# ------------- #
# model 1       #
# ------------- #

probs <- read_csv("Data/model1/probs.csv.gz")

lxs <-
  probs |> 
  filter(to > from) |> 
  mutate(transition = paste0(from,"-",to)) |> 
  select(-from, -to) |> 
  group_by(replicate, female, period5) |> 
  group_modify(~calc_occupancies(p_list=.x,
                                init=c(`1`=1,`2`=0),
                                delim="-",
                                age_interval = 0.25,
                                trans_col = "transition")) |> 
  rename(lx1 = `1`, lx2 = `2`) |> 
  mutate(lx1 = 4 * lx1,
         lx2 = 4 * lx2) 

write_csv(lxs, file = "Data/model1/lxs.csv.gz")

prev_stationary_summary <- 
  lxs |> 
  group_by(replicate,female, period5) |> 
  mutate(Lx1 = .25 * (lead(lx1,default = 0) + lx1) / 2,
         Lx2 = .25 * (lead(lx2,default = 0) + lx2) / 2,
         prev_st = Lx2 / (Lx1 + Lx2)) |> 
  group_by(female, period5, age) |> 
  summarize(prev_stationary = median(prev_st),
            lower = quantile(prev_st,0.025),
            upper = quantile(prev_st,0.975), .groups="drop")

write_csv(prev_stationary_summary,file = "Data/model1/prev_stationary_summary.csv.gz")

exs <-
  lxs |> 
  group_by(replicate,female, period5) |> 
  mutate(Lx1 = .25 * (lead(lx1,default = 0) + lx1) / 2,
         Lx2 = .25 * (lead(lx2,default = 0) + lx2) / 2) |> 
  summarize(DFLE = sum(Lx1),
            DLE = sum(Lx2), .groups="drop") 
  
write_csv(exs, file = "Data/model1/exs.csv.gz")

ex_summary <- 
  exs |> 
  ungroup() |> 
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, period5, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") 

write_csv(ex_summary,file = "Data/model1/exs_summary.csv.gz")

ex_summary |> 
  mutate(gender = if_else(female == 1,"women","men")) |> 
  ggplot() +
  geom_pointrange(aes(x = period5, y = e50, ymin = lower, ymax = upper, color = gender),position = position_dodge2(width=.3)) +
  facet_wrap(~state, scales = "free_y") +
  theme_minimal() +
  scale_color_manual(values = c(women = "#8b46ab", men = "#46ab6e"))

rm(lxs,exs,probs,ex_summary,prev_stationary_summary);gc()

# ------------- #
# model 2       #
# ------------- #

probs <- read_csv("Data/model2/probs.csv.gz")

lxs <-
  probs |> 
  filter(to > from) |> 
  mutate(transition = paste0(from,"-",to)) |> 
  select(-from, -to) |> 
  group_by(replicate, female, year) |> 
  group_modify(~calc_occupancies(p_list=.x,
                                 init=c(`1`=1,`2`=0),
                                 delim="-",
                                 age_interval = 0.25,
                                 trans_col = "transition")) |> 
  rename(lx1 = `1`, lx2 = `2`) |> 
  mutate(lx1 = 4 * lx1,
         lx2 = 4 * lx2) 

write_csv(lxs, file = "Data/model2/lxs.csv.gz")
prev_stationary_summary <- 
  lxs |> 
  group_by(replicate,female, year) |> 
  mutate(Lx1 = .25 * (lead(lx1,default = 0) + lx1) / 2,
         Lx2 = .25 * (lead(lx2,default = 0) + lx2) / 2,
         prev_st = Lx2 / (Lx1 + Lx2)) |> 
  group_by(female, year,age) |> 
  summarize(prev_stationary = median(prev_st),
            lower = quantile(prev_st,0.025),
            upper = quantile(prev_st,0.975), .groups="drop")

write_csv(prev_stationary_summary,file = "Data/model2/prev_stationary_summary.csv.gz")

exs <-
  lxs |> 
  group_by(replicate,female, year) |> 
  mutate(Lx1 = .25 * (lead(lx1,default = 0) + lx1) / 2,
         Lx2 = .25 * (lead(lx2,default = 0) + lx2) / 2) |> 
  summarize(DFLE = sum(Lx1),
            DLE = sum(Lx2), .groups="drop") 

write_csv(exs, file = "Data/model2/exs.csv.gz")

ex_summary <- 
  exs |> 
  ungroup() |> 
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, year, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") 
write_csv(ex_summary,file = "Data/model2/exs_summary.csv.gz")

ex_summary |> 
  mutate(gender = if_else(female == 1,"women","men")) |> 
  ggplot(aes(x=year, y = e50)) +
  geom_line(aes(color = gender)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = gender)) +
  facet_wrap(~state, scales = "free_y") +
  theme_minimal() +
  scale_color_manual(values = c(women = "#8b46ab", men = "#46ab6e")) +
  scale_fill_manual(values = c(women = "#8b46ab30", men = "#46ab6e30"))

rm(lxs,exs,probs,ex_summary,prev_stationary_summary);gc()


# ------------- #
# model 3       #
# ------------- #

probs <- read_csv("Data/model3/probs.csv.gz")

lxs <-
  probs |> 
  filter(to > from) |> 
  mutate(transition = paste0(from,"-",to)) |> 
  select(-from, -to) |> 
  group_by(replicate, female, year) |> 
  group_modify(~calc_occupancies(p_list=.x,
                                 init=c(`1`=1,`2`=0),
                                 delim="-",
                                 age_interval = 0.25,
                                 trans_col = "transition")) |> 
  rename(lx1 = `1`, lx2 = `2`) |> 
  mutate(lx1 = 4 * lx1,
         lx2 = 4 * lx2) 

write_csv(lxs, file = "Data/model3/lxs.csv.gz")

exs <-
  lxs |> 
  group_by(replicate,female, year) |> 
  mutate(Lx1 = .25 * (lead(lx1,default = 0) + lx1) / 2,
         Lx2 = .25 * (lead(lx2,default = 0) + lx2) / 2) |> 
  summarize(DFLE = sum(Lx1),
            DLE = sum(Lx2), .groups="drop") 

write_csv(exs, file = "Data/model3/exs.csv.gz")

prev_stationary_summary <- 
  lxs |> 
  group_by(replicate,female, year) |> 
  mutate(Lx1 = .25 * (lead(lx1,default = 0) + lx1) / 2,
         Lx2 = .25 * (lead(lx2,default = 0) + lx2) / 2,
         prev_st = Lx2 / (Lx1 + Lx2)) |> 
  group_by(female, year,age) |> 
  summarize(prev_stationary = median(prev_st),
            lower = quantile(prev_st,0.025),
            upper = quantile(prev_st,0.975), .groups="drop")
write_csv(prev_stationary_summary,file = "Data/model3/prev_stationary_summary.csv.gz")

ex_summary <- 
  exs |> 
  ungroup() |> 
  pivot_longer(c(DFLE,DLE), names_to = "state",values_to = "expectancy") |> 
  group_by(female, year, state) |> 
  summarize(e50 = median(expectancy),
            lower = quantile(expectancy,0.025),
            upper = quantile(expectancy,0.975), .groups = "drop") 
write_csv(ex_summary,file = "Data/model3/exs_summary.csv.gz")

  ex_summary |> 
  mutate(gender = if_else(female == 1,"women","men")) |> 
  ggplot(aes(x=year, y = e50)) +
  geom_line(aes(color = gender)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = gender)) +
  facet_wrap(~state, scales = "free_y") +
  theme_minimal() +
  scale_color_manual(values = c(women = "#8b46ab", men = "#46ab6e")) +
  scale_fill_manual(values = c(women = "#8b46ab30", men = "#46ab6e30"))


prev_stationary |> 
  mutate(gender = if_else(female == 1,"women","men")) |> 
  ggplot(aes(x=age, y = prev_stationary, color = year, group=year)) +
  geom_line(aes(color = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = year), color = "transparent", alpha = .2) +
  facet_wrap(~gender, scales = "free_y") +
  theme_minimal() +
  scale_color_binned_sequential() +
  scale_fill_binned_sequential()

  
rm(lxs,exs,probs,prev_stationary, ex_summary);gc()










