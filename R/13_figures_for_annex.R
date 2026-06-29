library(tidyverse)
library(grid)  # для стрелок
library(DiagrammeR)
library(ggplot2)
# Пример данных для 2 человек
df <- tibble(
  id = rep(1:2, each = 12),
  wave = rep(5:16, 2),
  observed = c(
    "No", "No", "Yes", "No", "Yes", "No", "Yes", "No", "No", "Yes", "No", "No",   # person 1
    "No", "No", "No", "Yes", "No", "No", "Yes", "No", "No", "No", "Yes", "No"      # person 2
  )
)

# Создаем absorbing state
df <- df %>%
  group_by(id) %>%
  mutate(absorbing = ifelse(observed == "Yes", 1, 0),
         absorbing = cummax(absorbing),
         absorbing_label = ifelse(absorbing == 1, "Yes", "No")) %>%
  ungroup()
# Высота линий для графика
df <- df %>%
  mutate(y_observed = id * 1.5,
         y_absorbing = y_observed - 0.5)

# График
ggplot(df, aes(x = wave)) +
  # Верхняя линия: фактические наблюдения
  geom_line(aes(y = y_observed, group = id), color = "black") +
  geom_point(aes(y = y_observed, color = observed), size = 4) +
  geom_text(aes(y = y_observed, label = observed), vjust = -1.2, size = 3.5) +
  
  # Нижняя линия: absorbing state
  geom_line(aes(y = y_absorbing, group = id), color = "grey30", linetype = "dashed") +
  geom_point(aes(y = y_absorbing, color = absorbing_label), size = 4) +
  geom_text(aes(y = y_absorbing, label = absorbing_label), vjust = -1.2, size = 3) +
  
  # Стрелка с пояснением
  annotate("segment", x = 5, xend = 16, y = 0.2, yend = 0.2, 
           arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
  annotate("text", x = 10.5, y = 0.4, 
           label = "All observations after first 'Yes' are treated as 'Yes'", 
           color = "red", size = 4, fontface = "bold") +
  
  scale_x_continuous(breaks = 5:16) +
  scale_y_continuous(breaks = c(1.5, 3), labels = c("Person 1", "Person 2")) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 11),
        legend.position = "top") +
  labs(x = "Wave", color = "Status")


ggsave(
  filename = "state_registration.jpeg",
  scale = 1,
  dpi = 300
)



grViz("
digraph dementia_states {
  
  # Общие стили графа
  graph [rankdir = LR, fontsize=12, labelloc='t', label='State-Space Diagram: Dementia and Death', fontname=Helvetica]
  
  # Узлы с цветами и формой
  node [shape = circle, style = filled, fontname = Helvetica, fontsize=12, width=1.2]
  
  NoDementia [label='No Dementia\\n(0)', fillcolor='lightblue']
  Dementia  [label='Dementia\\n(1)', fillcolor='orange']
  Death     [label='Death\\n(2)', fillcolor='red']
  
  # Переходы
  NoDementia -> Dementia [label='0→1', color='orange', fontcolor='orange', penwidth=2]
  NoDementia -> Death    [label='0→2', color='red', fontcolor='red', penwidth=2]
  Dementia -> Death      [label='1→2', color='red', fontcolor='red', penwidth=2]
  
  # Подпись снизу
  subgraph legend {
    label='Note: Dementia and Death are absorbing states'
    fontsize=10
    style=invis
  }
}
")

ggsave(
  filename = "state_space.jpeg",
  scale = 1,
  dpi = 300
)




haz_unadj <- read_csv("Data/model2/unadj_haz_replicates.csv.gz", show_col_types = FALSE) %>% 
  filter(year == 2019) %>% 
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
  unite(Transtion, c("from", "to"), sep = "-") %>% 
  mutate(Hazards = "Unadjusted")


haz <- read_csv("Data/model2/adj_haz_replicates.csv.gz", show_col_types = FALSE) %>%
  filter(year == 2019) %>%
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
  unite(Transtion, c("from", "to"), sep = "-") %>% 
  mutate(Hazards = "Adjusted")


# before and after adjustment hazards
haz %>% 
  full_join(haz_unadj) %>% 
  select(-year) %>% 
  mutate(female = ifelse(female == 1, "Female", "Male"),
         female = factor(female, levels = c("Male", "Female"))) %>% 
  separate(Transtion, c("From", "To")) %>% 
  ggplot(aes(x = age, y = p_med, color = To, linetype = Hazards)) +
  geom_line(linewidth = 1) + 
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = To),
    alpha = 0.12,
    color = NA
  ) +
  theme_bw(base_size = 14) +
  scale_y_log10()+
  facet_grid(female ~ From, switch = "y") + 
  theme(strip.placement  = "outside",
        strip.background = element_blank(),
        legend.position  = "bottom")

ggsave(
  filename = "state_registration.jpeg",
  scale = 1,
  dpi = 300
)








compare |>
  ggplot(aes(x = age, y = haz_med)) +
  # uncertainty bands (subtle!)
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = interaction(from, to)),
    alpha = 0.12,
    color = NA
  ) +
  # main hazard lines
  geom_line(
    aes(color = interaction(from, to),
        linetype = variant),
    linewidth = 0.9
  ) +
  # cleaner faceting (sex × period)
  facet_grid(female ~ period5,
             labeller = labeller(
               female = c(`0` = "Men", `1` = "Women")
             )) +
  # log scale
  scale_y_log10() +
  # colorblind-friendly palette
  # scale_color_brewer(palette = "Dark2", name = "Transition") +
  # scale_fill_brewer(palette = "Dark2", guide = "none") +
  # adjusted vs unadjusted
  scale_linetype_manual(
    values = c("solid", "dashed"),
    name = "Estimate",
    labels = c("Adjusted", "Unadjusted")
  ) +
  labs(
    x = "Age",
    y = "Hazard (log scale)",
    title = "Age-specific transition hazards",
    subtitle = "Adjusted vs. unadjusted estimates with 95% uncertainty intervals"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.spacing = unit(1.2, "lines")
  )

compare |>
  filter(period5 == "period 1") |>
  mutate(
    sex = if_else(female == 1, "Women", "Men"),
    transition = paste0(from, " \u2192 ", to),
    variant = factor(variant, levels = c("adjusted", "unadjusted"))
  ) |>
  ggplot(aes(x = age, y = haz_med)) +
  
  # subtle confidence bands
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = variant),
    alpha = 0.15,
    color = NA
  ) +
  
  # main lines
  geom_line(
    aes(color = variant, linetype = variant),
    linewidth = 1
  ) +
  
  # clean faceting: transitions × sex
  facet_wrap(sex ~ transition) +
  
  scale_y_log10() +
  
  # clean, publication-friendly colors
  # scale_color_manual(
  #   values = c("adjusted" = "#1b9e77", 
  #              "unadjusted" = "#d95f02"),
  #   name = "Estimate"
  # ) +
  # 
  # scale_fill_manual(
  #   values = c("adjusted" = "#1b9e77", 
  #              "unadjusted" = "#d95f02"),
  #   guide = "none"
  # ) +
  
  scale_linetype_manual(
    values = c("adjusted" = "solid", 
               "unadjusted" = "dashed"),
    guide = "none"
  ) +
  
  labs(
    x = "Age",
    y = "Hazard (log scale)",
    title = "Transition hazards by age (Period 1)",
    subtitle = "Comparison of adjusted and unadjusted estimates",
    caption = "Shaded areas indicate 95% uncertainty intervals"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "bottom",
    panel.spacing = unit(1.2, "lines")
  )









