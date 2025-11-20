
source("R/00_package_and_functions.R")

# load external mortality reference, used for adjustment
external_mort <- 
  read_csv("Data/external_ref_mc.csv.gz",show_col_types = FALSE)

# -----------------------------
# model 1
haz_unadj <- 
  read_csv("Data/model1/unadj_haz_replicates.csv.gz",show_col_types = FALSE)
prev <- 
  read_csv("Data/model1/prev_replicates.csv.gz",show_col_types = FALSE)

# model 1 has 5-year periods
external_mort_5 <-
  external_mort |> 
  mutate(period5 = case_when(
    between(year, 2004, 2010) ~ "period 1",
    between(year, 2010, 2015) ~ "period 2",
    year > 2015 ~ "period 3",
    TRUE ~ "other")) |> 
  filter(period5 != "other") |> 
  group_by(replicate,period5, female, age) |> 
  summarize(mx = mean(mx), .groups = "drop")


  
haz_adj_model1 <-
  haz_unadj |> 
  mutate(transition = paste0("h",from,"_",to)) |> 
  select(-from,-to) |> 
  pivot_wider(names_from = transition, values_from = haz) |> 
  left_join(prev, 
            by = join_by(replicate, age, female, period5)) |> 
  left_join(external_mort_5, 
            by = join_by(replicate, age, female, period5)) |> 
  mutate(
    Rx = h2_3 / h1_3,
    h1_3 = mx / (1 - prevalence + prevalence * Rx),
    h2_3 = h1_3 * Rx,
    h1_1 = -(h1_2 + h1_3),
    h2_2 = -(h2_1 + h2_3)) |> 
  select(replicate, age, female, period5, starts_with("h")) |> 
  pivot_longer(starts_with("h"), 
               names_to = "transition", 
               values_to = "haz") |> 
  mutate(transition = substr(transition,2,4)) |> 
  separate_wider_delim(transition,names=c("from","to"), delim="_") |> 
  mutate(from = as.integer(from),
         to = as.integer(to))
write_csv(haz_adj_model1,file = "Data/model1/adj_haz_replicates.csv.gz")

visualize_adjustment <- FALSE
if (visualize_adjustment){
unadj1 <-
haz_unadj |> 
  filter(to > from) |> 
  group_by(period5, female, age, from, to) |> 
  summarize(haz_med = median(haz),
            lower = quantile(haz,0.025),
            upper = quantile(haz, 0.975), .groups = "drop") |> 
  mutate(variant = "unadjusted")
adj1 <-
  haz_adj_model1 |> 
  filter(to > from) |> 
  group_by(period5, female, age, from, to) |> 
  summarize(haz_med = median(haz),
            lower = quantile(haz,0.025),
            upper = quantile(haz, 0.975), .groups = "drop") |> 
  mutate(variant = "adjusted") 

compare <- bind_rows(unadj1,
                     adj1)

compare |> 
  ggplot(aes(x =age, y = haz_med, color = interaction(from,to))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = interaction(from,to)),alpha=.3,color = "transparent") +
  facet_wrap(female~period5~variant) +
  scale_y_log10()
}

rm(external_mort_5)
rm(haz_unadj_model1)
rm(haz_unadj)
rm(prev)
rm(unadj1)
rm(haz_adj_model1)
gc()

# --------------------------------------
# model 2
# --------------------------------------

haz_unadj <- 
  read_csv("Data/model2/unadj_haz_replicates.csv.gz",show_col_types = FALSE)
prev <- 
  read_csv("Data/model2/prev_replicates.csv.gz",show_col_types = FALSE)

# pipeline too memory-greedy. Below this, find a low-memory data.table redux of 
# the same.

# haz_adj_model2 <-
#   haz_unadj |> 
#   mutate(transition = paste0("h",from,"_",to)) |> 
#   select(-from,-to) |> 
#   pivot_wider(names_from = transition, values_from = haz) |> 
#   left_join(prev, 
#             by = join_by(replicate, age, female, year)) |> 
#   left_join(external_mort, 
#             by = join_by(replicate, age, female, year)) |> 
#   mutate(
#     Rx = h2_3 / h1_3,
#     h1_3 = mx / (1 - prevalence + prevalence * Rx),
#     h2_3 = h1_3 * Rx,
#     h1_1 = -(h1_2 + h1_3),
#     h2_2 = -(h2_1 + h2_3)) |> 
#   select(replicate, age, female, year, starts_with("h")) |> 
#   pivot_longer(starts_with("h"), 
#                names_to = "transition", 
#                values_to = "haz") |> 
#   mutate(transition = substr(transition,2,4)) |> 
#   separate_wider_delim(transition,names=c("from","to"), delim="_") |> 
#   mutate(from = as.integer(from),
#          to = as.integer(to))


# Convert to data.table by reference (no copy on assignment)
setDT(haz_unadj)
setDT(prev)
setDT(external_mort)

# Keep only off-diagonals first (you don't need diag from the model)
haz_off <- haz_unadj[to > from]
rm(haz_unadj); gc()
# Create separate tables for each off-diagonal transition
h12 <- haz_off[from == 1 & to == 2,
               .(replicate, age, female, year, h1_2 = haz)]
h13 <- haz_off[from == 1 & to == 3,
               .(replicate, age, female, year, h1_3 = haz)]
h23 <- haz_off[from == 2 & to == 3,
               .(replicate, age, female, year, h2_3 = haz)]

# Drop haz_off and even haz_unadj to free memory
rm(haz_off); gc()


# Key everything on the join vars
setkey(h12, replicate, age, female, year)
setkey(h13, replicate, age, female, year)
setkey(h23, replicate, age, female, year)
setkey(prev, replicate, age, female, year)
setkey(external_mort, replicate, age, female, year)

# Build a compact group-level table: one row per (replicate, age, female, year)
haz_params <- h12[h13][h23]
rm(h12, h13, h23); gc()

# Add prevalence and mx IN PLACE (no massive copies)
haz_params[prev,         prevalence := i.prevalence]
haz_params[external_mort, mx        := i.mx]

# Now do your adjustment
haz_params[, Rx   := h2_3 / h1_3]
haz_params[, h1_3 := mx / (1 - prevalence + prevalence * Rx)]
haz_params[, h2_3 := h1_3 * Rx]
haz_params[, h1_1 := -(h1_2 + h1_3)]
haz_params[, h2_2 := -(h2_3)]
haz_params[, h2_1 := 0]
# Now re-expand to long hazard table (states 1 and 2)
haz_adj_model2 <- rbindlist(list(
  haz_params[, .(replicate, age, female, year, from = 1L, to = 1L, haz = h1_1)],
  haz_params[, .(replicate, age, female, year, from = 1L, to = 2L, haz = h1_2)],
  haz_params[, .(replicate, age, female, year, from = 1L, to = 3L, haz = h1_3)],
  haz_params[, .(replicate, age, female, year, from = 2L, to = 2L, haz = h2_1)],
  haz_params[, .(replicate, age, female, year, from = 2L, to = 2L, haz = h2_2)],
  haz_params[, .(replicate, age, female, year, from = 2L, to = 3L, haz = h2_3)]
), use.names = TRUE)

rm(haz_params);rm(prev);gc()
write_csv(haz_adj_model2,file = "Data/model2/adj_haz_replicates.csv.gz")
rm(haz_adj_model2);gc()

# --------------------------------------
# model 3
# --------------------------------------

haz_unadj <- 
  read_csv("Data/model3/unadj_haz_replicates.csv.gz",show_col_types = FALSE)
prev <- 
  read_csv("Data/model3/prev_replicates.csv.gz",show_col_types = FALSE)


setDT(haz_unadj)
setDT(prev)
setDT(external_mort)

# Keep only off-diagonals first (you don't need diag from the model)
haz_off <- haz_unadj[to > from]
rm(haz_unadj); gc()
# Create separate tables for each off-diagonal transition
h12 <- haz_off[from == 1 & to == 2,
               .(replicate, age, female, year, h1_2 = haz)]
h13 <- haz_off[from == 1 & to == 3,
               .(replicate, age, female, year, h1_3 = haz)]
h23 <- haz_off[from == 2 & to == 3,
               .(replicate, age, female, year, h2_3 = haz)]

# Drop haz_off and even haz_unadj to free memory
rm(haz_off); gc()


# Key everything on the join vars
setkey(h12, replicate, age, female, year)
setkey(h13, replicate, age, female, year)
setkey(h23, replicate, age, female, year)
setkey(prev, replicate, age, female, year)
setkey(external_mort, replicate, age, female, year)

# Build a compact group-level table: one row per (replicate, age, female, year)
haz_params <- h12[h13][h23]
rm(h12, h13, h23); gc()

# Add prevalence and mx IN PLACE (no massive copies)
haz_params[prev,         prevalence := i.prevalence]
haz_params[external_mort, mx        := i.mx]

# Now do your adjustment
haz_params[, Rx   := h2_3 / h1_3]
haz_params[, h1_3 := mx / (1 - prevalence + prevalence * Rx)]
haz_params[, h2_3 := h1_3 * Rx]
haz_params[, h1_1 := -(h1_2 + h1_3)]
haz_params[, h2_2 := -(h2_3)]
haz_params[, h2_1 := 0]
# Now re-expand to long hazard table (states 1 and 2)
haz_adj_model3 <- rbindlist(list(
  haz_params[, .(replicate, age, female, year, from = 1L, to = 1L, haz = h1_1)],
  haz_params[, .(replicate, age, female, year, from = 1L, to = 2L, haz = h1_2)],
  haz_params[, .(replicate, age, female, year, from = 1L, to = 3L, haz = h1_3)],
  haz_params[, .(replicate, age, female, year, from = 2L, to = 2L, haz = h2_1)],
  haz_params[, .(replicate, age, female, year, from = 2L, to = 2L, haz = h2_2)],
  haz_params[, .(replicate, age, female, year, from = 2L, to = 3L, haz = h2_3)]
), use.names = TRUE)

rm(haz_params);rm(prev);gc()
write_csv(haz_adj_model3,file = "Data/model3/adj_haz_replicates.csv.gz")
rm(haz_adj_model3);gc()
