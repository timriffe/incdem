# change colors to purple and green

library(tidyverse)
library(scales)
haz  <- read_csv("Data/model2/adj_haz_replicates.csv.gz")
prev <- read_csv("Data/model2/prev_replicates.csv.gz")


prev |> 
  filter(year %in% c(2004, 2019)) %>% 
  group_by(female, year, age) |> 
  summarize(prev_median = median(prevalence),
            lower       = quantile(prevalence, .025),
            upper       = quantile(prevalence, .975), 
            .groups = "drop") |>
  mutate(year = as.factor(year),
         female = factor(
           female,
           levels = c(0, 1),
           labels = c("Male", "Female"))) %>% 
  filter(age == 50)


# start with prevalence
prev |>
  filter(year %in% c(2004, 2019)) %>%
  group_by(female, year, age) |>
  summarize(prev_median = median(prevalence),
            lower       = quantile(prevalence, .025),
            upper       = quantile(prevalence, .975), 
            .groups = "drop") |>
  mutate(year = as.factor(year),
           female = factor(
             female,
             levels = c(0, 1),
             labels = c("Male", "Female"))) %>% 
  ggplot(aes(x = age, y = prev_median, color = year)) +
  geom_line(linewidth = 0.9) +
  geom_ribbon(mapping = aes(ymin = lower, 
                            ymax = upper, 
                            fill = year),
              alpha = .3,
              linewidth = 0.1,
              linetype = "solid")+
  facet_grid( ~ female, switch = "y") + 
  scale_x_continuous(
    breaks = seq(50, 100, 5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_brewer(
    palette = "Dark2"
  )+
  scale_color_brewer(
    palette = "Dark2"
  ) +
  labs(
    x = "Age",
    y = "Prevalence of dementia",
    color = "Year",
    fill = "Year"
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  coord_cartesian(clip = "on") +
  theme_bw(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.1, "lines"),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    strip.text = element_text(
      face = "bold",
      size = 11, 
      color = "black"
    ),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box = "horizontal"
  )

ggsave(filename = "figures_paper/fig1_prevalence.jpeg", scale = 1)





fig2 <- haz %>% 
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


fig2 |>
  filter(year %in% c(2004, 2019)) %>% 
  mutate(
    female = factor(
      female,
      levels = c(0, 1),
      labels = c("Male", "Female")
    ),
    year = factor(year),
    Transtion = case_when(
      Transtion == "H-D" ~ "Dementia-free to Death",
      Transtion == "U-D" ~ "Dementia to Death",
      Transtion == "H-U" ~ "Dementia-free to Dementia"
    ),
    Transtion = factor(Transtion, levels = c("Dementia-free to Dementia",
                                             "Dementia-free to Death",
                                             "Dementia to Death")
                       )) |>
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
    linewidth = 0.1,
    linetype = "solid",
    show.legend = FALSE
  ) +
  scale_fill_brewer(
    palette = "Dark2"
  )+
  scale_color_brewer(
    palette = "Dark2"
  ) +
  geom_line(linewidth = 0.9) +
  
  facet_grid(
    female ~ Transtion,
    switch = "y"
  ) +
  
  scale_y_log10(
    # labels = scales::label_scientific(digits = 1)
  ) +
  
  scale_x_continuous(
    breaks = seq(50, 100, 10),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    x = "Age",
    y = "Transition hazard (log scale)",
    color = "Year",
    fill = "Year"
  ) +
  
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  
  coord_cartesian(clip = "on") +
  
  theme_bw(base_size = 12) +
  
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.1, "lines"),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    strip.text = element_text(
      face = "bold",
      size = 10, 
      color = "black"
    ),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", color = "black"),
    legend.box = "horizontal",
    plot.title = element_text(
      face = "bold",
      size = 14, 
      color = "black"
    ))

ggsave(filename = "figures_paper/fig2_hazards.jpeg", scale = 1.2)
# ------------------------------------------------------------------- #
# probs2 <- read_csv("Data/model2/probs.csv.gz")
# 
# probs2 |> 
#   filter(to > from) |> 
#   group_by(period5, female, age, from, to) |> 
#   summarize(p_med = median(p),
#             lower = quantile(p, 0.025),
#             upper = quantile(p, 0.975)) |> 
#   ggplot(aes(x = age, y = p_med, color = interaction(from, to))) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = interaction(from, to)), color = "transparent", alpha = .3) +
#   facet_grid(female ~ period5, switch = "y") +
#   scale_y_log10() + 
#   theme(strip.placement  = "outside",
#         strip.background = element_blank())



# ------------------------------------------------------------------- #
# LE with and without dimentia

ex2 <- read_csv("Data/model2/e50.csv.gz")


ex2 %>% 
  filter(year %in% c(2004, 2019))


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
    linewidth = 0.1,
    linetype = "solid"
  ) +
  # central estimate
  geom_line(linewidth = 1) +
  scale_fill_brewer(
    palette = "Dark2"
  )+
  scale_color_brewer(
    palette = "Dark2"
  ) +
  geom_point(size = 1.8) +
  facet_wrap(
    State ~ .,
    # switch = "y",
    scales = "free_y",
  )+
  scale_x_continuous(
    breaks = seq(2004, 2019, 3)
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
    # strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11, 
      color = "black"
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", color = "black"),
    panel.spacing = unit(1.2, "lines"))

ggsave(filename = "figures_paper/fig3_e50.jpeg", scale = 1)

# ------------------------------------------------------------------- #
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
  filter(year %in% c(2004, 2019))


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
    linewidth = 0.1,
    linetype = "solid"
  ) +
  scale_fill_brewer(
    palette = "Dark2"
  )+
  scale_color_brewer(
    palette = "Dark2"
  ) +
  # central estimate
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  facet_grid(
    ~ female
  ) +
  scale_x_continuous(
    breaks = seq(2004, 2019, 3)
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
      size = 11, 
      color = "black"
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", color = "black"),
    panel.spacing = unit(1.2, "lines")
  )

ggsave(filename = "figures_paper/fig4_prop.jpeg", scale = 1)
