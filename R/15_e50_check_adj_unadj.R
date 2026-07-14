# males line 605
# females 295
source("R/00_dependencies.R")
source("R/01_functions.R")
source("R/02_prepare_hrs.R")
source("R/01_functions.R")
source("R/01_functions_extra.R")

hrs_to_fit_prepped <- 
  hrs_to_fit |>
  mutate( year = floor(obs_date)) |>
  select(hhidpn, female, education,
         pwt = wtcrnh, age, hypertension:stroke,
         state_msm, ever_dementia, obstype,
         year) |> 
  filter(year > 2003)

Q_mod2 <- rbind(
  c(-0.2,  0.1, 0.1),
  c( 0,   -0.01, 0.1),
  c( 0,    0,   0)
)

haz_unadj <- hrs_to_fit_prepped |>
  fit_msm(
    strat_vars    = c("female"),
    covariate_var = c("age", "year"),
    age_from_to   = c(50, 100),
    age_int       = 0.25,
    spline_df     = NULL,
    spline_type   = "ns",
    Q             = Q_mod2
  ) |>
  dplyr::mutate(replicate = 1L, .before = 1)

prev <- hrs_to_fit_prepped |>
  fit_prev(
    strat_vars      = c("female"),
    covariate_var   = c("age", "year"),
    age_from_to     = c(50, 100),
    age_int         = 0.25,
    condition_state = 2L,
    spline_df       = NULL,
    spline_type     = "ns",
    state_var       = "state_msm",
    weight_col      = "pwt",
    id_col          = "hhidpn",
    exclude_state   = 3L
  ) |>
  dplyr::mutate(replicate = 1L, .before = 1)

plan(multisession, workers = 10)

probs_unadj <-
  haz_unadj |> 
  hazards_to_discrete(
    age_interval = 0.25,
    id_cols      = c("age", "female", "year"),
    n_cores      = 10,
    parallel     = "future")

e50_unadj <-
  probs_unadj |> 
  mutate(transition = paste(from,to,sep="-")) |> 
  calc_exs(from_age = 50, 
           age_interval = 0.25, 
           init = c(`1` = 1, `2` = 0) ,
           init_method = "init", 
           from_col = "from", 
           to_col = "to", 
           age_col = "age", 
           p_col = "p", 
           group_cols = c("female", "year"), 
           trans_col = "transition") |> 
  mutate(LE = DFLE + DLE,
         adjusted = FALSE)

e50_unadj_lt <-
  haz_unadj |> 
  filter(to == 3,
         from < 3) |> 
  mutate(transition = if_else(from == 1,"HD","UD")) |> 
  select(-from,-to,-replicate) |> 
  pivot_wider(names_from = transition, values_from = haz) |> 
  left_join(prev |> select(-replicate),
            by = join_by(female, age, year)) |> 
  mutate(mx = HD * (1 - prevalence) + UD * prevalence) |> 
  group_by(female, year) |> 
  mutate(ax = .125,
         qx = (.25 * mx) / (1 + (.25 - ax) * mx),
         lx = c(1,cumprod(1-qx[-n()])),
         dx = lx * qx,
         Lx = lx * .25 - (1 - ax) * dx) |> 
  summarize(LE = sum(Lx))

external_mort <- 
  read_csv("Data/external_ref_point.csv",show_col_types = FALSE)

haz_adj <-
  haz_unadj |> 
  mutate(female = as.integer(female)-1) |> 
  mutate(transition = paste0("h",from,"_",to)) |> 
  select(-from,-to) |> 
  pivot_wider(names_from = transition, values_from = haz) |> 
  left_join(prev |> select(-replicate) |> mutate(female = as.integer(female)-1), 
            by = join_by(age, female, year)) |> 
  left_join(external_mort, 
            by = join_by(age, female, year)) |> 
  mutate(
    Rx = h2_3 / h1_3,
    h1_3 = mx / (1 - prevalence + prevalence * Rx),
    h2_3 = h1_3 * Rx,
    h1_1 = -(h1_2 + h1_3),
    h2_2 = -(h2_1 + h2_3)) |> 
  select(replicate, age, female, year, starts_with("h")) |> 
  pivot_longer(starts_with("h"), 
               names_to = "transition", 
               values_to = "haz") |> 
  mutate(transition = substr(transition,2,4)) |> 
  separate_wider_delim(transition,names=c("from","to"), delim="_") |> 
  mutate(from = as.integer(from),
         to = as.integer(to)) |> 
  filter(from < 3)

haz_compare <-
  haz_adj |> 
  mutate(adjusted = TRUE) |> 
  bind_rows(haz_unadj |> mutate(female = as.integer(female)-1,adjusted = FALSE))

# in 2005, the adjustment increases HD and UD, but in 2019 it decreases them
# I didn't expect HRS to overestimate mortality ever..
haz_compare |> 
  filter(to > from, year %in% c(2005,2019)) |> 
  ggplot(aes(x = age, y = haz, color = interaction(from, to), linetype = adjusted)) +
  geom_line() + 
  facet_wrap(year~female) +
  scale_y_log10() +
  geom_line(data = external_mort |> filter(year %in% c(2005,2019)) |> rename(haz=mx) |> mutate(adjusted = FALSE), color = "black")

plan(multisession, workers = 10)
probs_adj <-
  haz_adj |> 
  arrange(year,female,age, from,to) |> 
  hazards_to_discrete(
    age_interval = 0.25,
    id_cols      = c("age", "female", "year"),
    n_cores      = 10,
    parallel     = "future")

e50_adj <-
  probs_adj |> 
  mutate(transition = paste(from,to,sep="-")) |> 
  calc_exs(from_age = 50, 
           age_interval = 0.25, 
           init = c(`1` = 1, `2` = 0) ,
           init_method = "init", 
           from_col = "from", 
           to_col = "to", 
           age_col = "age", 
           p_col = "p", 
           group_cols = c("female", "year"), 
           trans_col = "transition") |> 
  mutate(LE = DFLE + DLE,
         adjusted = TRUE)

e50_adj_lt <-
  haz_adj |> 
  filter(to == 3,
         from < 3) |> 
  mutate(transition = if_else(from == 1,"HD","UD")) |> 
  select(-from,-to,-replicate) |> 
  pivot_wider(names_from = transition, values_from = haz) |> 
  left_join(prev |> select(-replicate) |> mutate(female = as.integer(female) - 1),
            by = join_by(female, age, year)) |> 
  mutate(mx = HD * (1 - prevalence) + UD * prevalence) |> 
  group_by(female, year) |> 
  mutate(ax = .125,
         qx = (.25 * mx) / (1 + (.25 - ax) * mx),
         lx = c(1,cumprod(1-qx[-n()])),
         dx = lx * qx,
         Lx = lx * .25 - (1 - ax) * dx) |> 
  summarize(LE = sum(Lx), .groups= "drop") |> 
  mutate(method = "lt",
         adjusted = TRUE)

e50 <-
  bind_rows(e50_adj|> 
              mutate(method = "mslt",
                     adjusted = TRUE), 
            e50_unadj |>  
              mutate(method = "mslt",
                     adjusted = FALSE,
                     female = as.integer(female)-1),
            e50_adj_lt|> 
              mutate(method = "lt",
                     adjusted = TRUE), 
            e50_unadj_lt |>  
              mutate(method = "lt",
                     adjusted = FALSE,
                     female = as.integer(female)-1)) 
  

# It seems to be the case that 
e50 |> 
  ggplot(aes(x= year, 
             y = LE, 
             color = adjusted, 
             linetype = method)) +
  facet_wrap(~female) +
  geom_point() + 
  geom_line()
