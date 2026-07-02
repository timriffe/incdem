# males line 605
# females 295
# source("R/00_dependencies.R")
# source("R/01_functions_extra.R")
# # which replicate is the median e50?
# probs <- read_csv("Data/model2/probs.csv.gz")
# # # expectancies, age 50
# e50 <- calc_exs(probs = probs,
#                 from_age = 50,
#                 age_interval = 0.25,
#                 init = c(`1` = 1, `2` = 0),
#                 init_method = "init",
#                 from_col = "from",
#                 to_col = "to",
#                 age_col = "age",
#                 p_col = "p",
#                 group_cols = c("replicate", "female", "year"),
#                 trans_col = "transition") |>
#   ungroup()


# females 295
# males 605
df <- e50 %>%
  group_by(year, female) %>% 
  # filter(year == 2019) %>% 
  mutate(LE = DFLE + DLE,
         res = LE - median(LE)) %>% 
  filter(res == min(abs(res))) %>% 
  select(1, 2, 3)

# real tibble of replicates and corresponding values
df <- tribble(
  ~replicate, ~female, ~year,
  155, 1, 2010,
  156, 0, 2006,
  183, 0, 2007,
  208, 1, 2004,
  256, 0, 2015,
  283, 1, 2005,
  295, 1, 2019,
  344, 0, 2010,
  369, 1, 2013,
  381, 0, 2011,
  415, 1, 2015,
  550, 1, 2011,
  596, 1, 2006,
  605, 0, 2019,
  618, 0, 2009,
  680, 1, 2017,
  753, 0, 2014,
  767, 0, 2017,
  768, 0, 2004,
  779, 0, 2012,
  815, 1, 2008,
  826, 1, 2014,
  869, 1, 2007,
  903, 1, 2009,
  960, 0, 2013,
  982, 0, 2018
)


prev          <- read_csv("Data/model2/prev_replicates.csv.gz",      show_col_types = FALSE)
hazards_unadj <- read_csv("Data/model2/unadj_haz_replicates.csv.gz", show_col_types = FALSE)
hazards_adj   <- read_csv("Data/model2/adj_haz_replicates.csv.gz",   show_col_types = FALSE)

prev1 <- prev %>% 
  semi_join(
    df,
    by = c("replicate", "female", "year")
  )  
  # filter(female == 0 & replicate == 605 |
  #        female == 1 & replicate == 295,
  #        year == 2019) %>% 
  # select(-year)

unadj <- hazards_unadj %>% 
  semi_join(
    df,
    by = c("replicate", "female", "year")
  )  %>% 
  filter(
    # female == 0 & replicate == 605 |
    #   female == 1 & replicate == 295,
    #      year == 2019,
         to > from) |> # create c from and to columns
  mutate(from = ifelse(from == 1, "H", "U"),
         to   = case_when(
           to == 1 ~ "H",
           to == 2 ~ "U",
           to == 3 ~ "D"
         )) %>% 
  mutate(transition = paste0(from, "_", to)) %>%
  select(-from, -to) %>%
  pivot_wider(names_from  = transition, 
              values_from = haz) %>%
  select(-H_U) %>%
  left_join(
    prev,
    by = join_by(replicate, age, female, year)
  ) %>% 
  mutate(
    mx = (1 - prevalence) * H_D + prevalence  * U_D,
    female = ifelse(female == 1, "Women", "Men")) %>%
  arrange(replicate, female, year, age) %>% 
  select(-replicate) %>% 
  mutate(Hazards = "Before adjustment")%>% 
  select(-c(H_D, U_D, prevalence))

# ABSOLUTELY SAME CODE
adj <- hazards_adj %>% 
  semi_join(
    df,
    by = c("replicate", "female", "year")
  )  %>%
  filter(
    # female == 0 & replicate == 605 |
    #   female == 1 & replicate == 295,
    #      year == 2019,
    to > from) |> # create c from and to columns
  mutate(from = ifelse(from == 1, "H", "U"),
         to   = case_when(
           to == 1 ~ "H",
           to == 2 ~ "U",
           to == 3 ~ "D"
         )) %>% 
  mutate(transition = paste0(from, "_", to)) %>%
  select(-from, -to) %>%
  pivot_wider(names_from  = transition, 
              values_from = haz) %>%
  select(-H_U) %>%
  left_join(
    prev,
    by = join_by(replicate, age, female, year)
  ) %>% 
  mutate(
    mx = (1 - prevalence) * H_D + prevalence  * U_D,
    female = ifelse(female == 1, "Women", "Men")) %>%
  arrange(replicate, female, year, age) %>% 
  select(-replicate) %>% 
  mutate(Hazards = "After adjustment")%>% 
  select(-c(H_D, U_D, prevalence))

# ------------------------------------------------------------------- #
# e50
z <- adj %>% 
  full_join(unadj) %>%
  group_by(Hazards, female, year) %>% 
  summarise(
    # I assume quarters     
    e50 = sum(exp(-cumsum(mx * 0.25))) * 0.25,
    .groups = "drop"
  )

# Something wrong
z %>% 
  ggplot(
  aes(
    x = year,
    y = e50,
    color = female,
    linetype = Hazards,
    group = interaction(female, Hazards)
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 2
  )
