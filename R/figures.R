library(tidyverse)



hrs_to_fit_prepped$year %>% unique() %>% sort()


# т.е. мы с 7-й волны берем так как раньше мало наблюдений просто
# from 2004
# В самом начале второго файла он с 2004 года фильтрует я не ебу почему, но вот так
# hrs_to_fit_prepped вот тут вот прям в начале
# надо поправить
haz <- read_csv("Data/model3/adj_haz_replicates.csv.gz")

fig1 <- haz %>% 
  filter(to > from) |> # create c from and to columns
  summarize(p_med = median(haz),
            lower = quantile(haz, 0.025),
            upper = quantile(haz, 0.975), 
            .by = c(year, female, age, from, to)) %>% 
  mutate(from = ifelse(from == 1, "H", "U"),
         to   = case_when(
           to == 1 ~ "H",
           to == 2 ~ "U",
           to == 3 ~ "D"
         )) %>% 
  unite(Transtion, c("from", "to"), sep = "-")


fig1 %>% 
  select(-c(lower, upper)) %>% 
  pivot_wider(names_from  = Transtion,
              values_from = p_med) %>% 
  mutate(Ra = `U-D` / `H-D`) %>% 
  select(-c(`H-U`, `H-D`, `U-D`))


# z |>
#   ggplot(aes(x = age, y = p_med, color = Transtion)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = Transtion), 
#               color = "transparent", 
#               alpha = .3) +
#   facet_grid(female ~ year, switch = "y") +
#   scale_y_log10() + 
#   theme(strip.placement  = "outside",
#         strip.background = element_blank(),
#         legend.position = "bottom")



# fig1 |>
#   filter(year %in% c(2004, 2009, 2014, 2019)) %>% 
#   mutate(
#     female = factor(female,
#                     levels = c(0, 1),
#                     labels = c("Male", "Female")),
#     year = factor(year),
#     Transtion = factor(
#       Transtion,
#       levels = c("H-U", "U-H", "H-D", "U-D")
#     )
#   ) |>
#   ggplot(aes(x = age,
#              y = p_med,
#              color = Transtion,
#              fill = Transtion)) +
#   # confidence bands
#   geom_ribbon(aes(ymin = lower,
#                   ymax = upper),
#               alpha = 0.18,
#               linewidth = 0,
#               show.legend = FALSE) +
#   # central estimate
#   geom_line(linewidth = 0.9) +
#   facet_grid(
#     female ~ year,
#     switch = "y"
#   ) +
#   scale_y_log10(
#     labels = scales::label_scientific(digits = 1)
#   ) +
#   scale_x_continuous(
#     breaks = seq(50, 100, 10),
#     expand = expansion(mult = c(0.01, 0.02))
#   ) +
#   labs(
#     x = "Age",
#     y = "Transition hazard (log scale)",
#     color = "Transition"
#     ) +
#   guides(
#     color = guide_legend(
#       nrow = 1,
#       byrow = TRUE
#     )
#   ) +
#   coord_cartesian(clip = "off") +
#   theme_bw(base_size = 12) +
#   theme(
#     strip.placement = "outside",
#     strip.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.spacing = unit(1.1, "lines"),
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 10),
#     strip.text = element_text(
#       face = "bold",
#       size = 11
#     ),
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold"),
#     legend.box = "horizontal",
#     plot.title = element_text(
#       face = "bold",
#       size = 14
#     ),
#     plot.subtitle = element_text(size = 11)
#   )
# 
# ggsave(filename = "fig1.jpeg", scale = 1)

fig1 |>
  filter(year %in% c(2004, 2019)) %>% 
  mutate(
    female = factor(
      female,
      levels = c(0, 1),
      labels = c("Male", "Female")
    ),
    year = factor(year),
    Transtion = factor(
      Transtion,
      levels = c("H-U", "U-H", "H-D", "U-D")
    )
  ) |>
  ggplot(
    aes(
      x = age,
      y = p_med,
      color = year,
      fill = year
    )
  ) +
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.15,
    linewidth = 0,
    show.legend = FALSE
  ) +
  
  geom_line(linewidth = 0.9) +
  
  facet_grid(
    female ~ Transtion,
    switch = "y"
  ) +
  
  scale_y_log10(
    labels = scales::label_scientific(digits = 1)
  ) +
  
  scale_x_continuous(
    breaks = seq(50, 100, 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  
  labs(
    x = "Age",
    y = "Transition hazard (log scale)",
    color = "Year"
  ) +
  
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 12) +
  
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.1, "lines"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal",
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(size = 11)
  )

ggsave(filename = "fig1_1.jpeg", scale = 1)


# haz1 <- read_csv("Data/model1/adj_haz_replicates.csv.gz")
# haz2 <- read_csv("Data/model2/adj_haz_replicates.csv.gz")
# haz3 <- read_csv("Data/model3/adj_haz_replicates.csv.gz")


# probs1 <- read_csv("Data/model1/probs.csv.gz")
# probs2 <- read_csv("Data/model2/probs.csv.gz")
probs2 <- read_csv("Data/model2/probs.csv.gz")

# ------------------------------------------------------------------- #

probs1
probs2
probs3


probs2 |> 
  filter(to > from) |> 
  group_by(period5, female, age, from, to) |> 
  summarize(p_med = median(p),
            lower = quantile(p, 0.025),
            upper = quantile(p, 0.975)) |> 
  ggplot(aes(x = age, y = p_med, color = interaction(from, to))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = interaction(from, to)), color = "transparent", alpha = .3) +
  facet_grid(female ~ period5, switch = "y") +
  scale_y_log10() + 
  theme(strip.placement  = "outside",
        strip.background = element_blank())



# ------------------------------------------------------------------- #
# LE with and without dimentia
# we need tyo adjust the 9 expectancies code for figure 2, but should not be difficult

ex2 <- read_csv("Data/model2/e50.csv.gz")

ex2 |>
  mutate(
    state = ifelse(
      state == "DFLE",
      "Dementia free",
      "With dementia"
    ),
    female = factor(
      female,
      levels = c(0, 1),
      labels = c("Men", "Women")
    ),
    State = factor(
      state,
      levels = c("Dementia free", "With dementia")
    )
  ) |>
  ggplot(
    aes(
      x = year,
      y = e50,
      color = female,
      fill = female
    )
  ) +
  # uncertainty intervals
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.18,
    linewidth = 0
  ) +
  # central estimate
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_wrap(
    ~ State,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = seq(2004, 2019, 2)
  ) +
  labs(
    x = "Year",
    y = "Life expectancy at age 50",
    color = "Sex",
    fill = "Sex"
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.spacing = unit(1.2, "lines")
  )

ggsave(filename = "fig2.jpeg", scale = 1)


# ------------------------------------------------------------------- #
# we will need to slightly update the conditional proportions figure
# reshape
# prop <- ex2  %>%
#   pivot_wider(
#     names_from  = state,
#     values_from = c(e50, lower, upper)
#   ) %>%
#   set_names(str_remove(names(.), "e50_")) %>%
#   mutate(
#
#     # central estimate
#     prop = DLE / (DLE + DFLE),
#
#     # approximate CI
#     # lower dementia share:
#     # smallest possible DLE / largest possible total
#     prop_lower =
#       lower_DLE / (upper_DFLE + lower_DLE),
#
#     # upper dementia share:
#     # largest possible DLE / smallest possible total
#     prop_upper =
#       upper_DLE / (lower_DFLE + upper_DLE),
#
#     female = factor(
#       female,
#       levels = c(0, 1),
#       labels = c("Men", "Women")
#     )
#   )



# # figure
# ggplot(
#   prop,
#   aes(
#     x = year,
#     y = prop,
#     color = female,
#     fill = female
#   )
# ) +
#
#   geom_ribbon(
#     aes(
#       ymin = prop_lower,
#       ymax = prop_upper
#     ),
#     alpha = 0.18,
#     linewidth = 0
#   ) +
#
#   geom_line(linewidth = 1) +
#
#   geom_point(size = 2) +
#
#   scale_y_continuous(
#     labels = scales::percent_format(accuracy = 1)
#   ) +
#
#   scale_x_continuous(
#     breaks = seq(2004, 2019, 2)
#   ) +
#
#   labs(
#     x = "Year",
#     y = "Proportion of remaining life lived with dementia",
#     color = "Sex",
#     fill = "Sex"
#   ) +
#
#   guides(
#     color = guide_legend(
#       nrow = 1,
#       byrow = TRUE
#     )
#   ) +
#
#   theme_bw(base_size = 12) +
#
#   theme(
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 10),
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold")
#   )
# 
# ggsave(filename = "fig3.jpeg", scale = 1)



prop <- read_csv("Data/model2/prop50.csv.gz")

prop |>
  mutate(
    female = factor(
      female,
      levels = c(0, 1),
      labels = c("Men", "Women")
    )
  ) |>
  mutate(across(c(prop50, lower, upper), ~ 1 - .)) %>% 
  ggplot(
    aes(
      x = year,
      y = prop50,
      color = female,
      fill = female
    )
  ) +
  # uncertainty intervals
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.18,
    linewidth = 0
  ) +
  # central estimate
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_wrap(
    ~ female,
    ncol = 1
  ) +
  scale_x_continuous(
    breaks = seq(2004, 2019, 2)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Year",
    y = "Proportion of remaining life spent with dementia at age 50",
    color = "Sex",
    fill = "Sex"
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),

    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.spacing = unit(1.2, "lines")
  )

ggsave(filename = "fig3.jpeg", scale = 1)

