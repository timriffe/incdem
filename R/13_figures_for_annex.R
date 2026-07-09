library(tidyverse)
library(grid) # arrows
library(DiagrammeR)
library(ggplot2)

# ------------------------------------------------------------------------ #
# example for 2 people
df <- tibble(
  id = rep(1:2, each = 12),
  wave = rep(5:16, 2),
  observed = c(
    "No", "No", "Yes", "No", "Yes", "No", "Yes", "No", "No", "Yes", "No", "No",   # person 1
    "No", "No", "No", "Yes", "No", "No", "Yes", "No", "No", "No", "Yes", "No"      # person 2
  )
)

# absorbing state
df <- df %>%
  group_by(id) %>%
  mutate(absorbing = ifelse(observed == "Yes", 1, 0),
         absorbing = cummax(absorbing),
         absorbing_label = ifelse(absorbing == 1, "Yes", "No")) %>%
  ungroup()

# line height
df <- df %>%
  mutate(y_observed = id * 1.5,
         y_absorbing = y_observed - 0.5)

# plot
ggplot(df, aes(x = wave)) +
  # upper line actual obs
  geom_line(aes(y = y_observed, group = id), color = "black") +
  geom_point(aes(y = y_observed, color = observed), size = 4) +
  geom_text(aes(y = y_observed, label = observed), vjust = -1.2, size = 3.5) +
  
  # lower line absorb state
  geom_line(aes(y = y_absorbing, group = id), color = "grey30", linetype = "dashed") +
  geom_point(aes(y = y_absorbing, color = absorbing_label), size = 4) +
  geom_text(aes(y = y_absorbing, label = absorbing_label), vjust = -1.2, size = 3) +
  
  # arrow with explaination
  annotate("segment", x = 5, xend = 16, y = 0.2, yend = 0.2, 
           arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
  annotate("text", x = 10.5, y = 0.4, 
           label = "All observations after first 'Yes (Dementia)' are treated as 'Yes'", 
           color = "red", size = 4, fontface = "bold") +
  scale_x_continuous(breaks = 5:16) +
  scale_y_continuous(breaks = c(1.5, 3), labels = c("Person 1", "Person 2")) +
  theme_minimal(base_family = "Arial", base_size = 10) +
  theme(  axis.title = element_text(size = 10, color = "black", face = "bold"),
          axis.text = element_text(size = 10, color = "black"),
          strip.text = element_text(size = 10, face = "bold", color = "black"),
          legend.title = element_text(size = 10, face = "bold", color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          legend.position = "bottom") +
  labs(x     = "HRS Wave",
       color = "Dementia status",
       y     = "Observed person") + 
  scale_fill_brewer(
    palette = "Dark2"
  )+
  scale_color_brewer(
    palette = "Dark2"
  )

ggsave(
  filename = "figures_annex/state_registration.pdf",
  device = cairo_pdf,
  dpi = 300
)

# ------------------------------------------------------------------------ #
# State space
# grViz("
# digraph dementia_states {
#   
#   # Общие стили графа
#   graph [rankdir = LR, fontsize=12, labelloc='t', label='State-Space Diagram: Dementia and Death', fontname=Helvetica]
#   
#   # Узлы с цветами и формой
#   node [shape = circle, style = filled, fontname = Helvetica, fontsize=12, width=1.2]
#   
#   NoDementia [label='No Dementia\\n(0)', fillcolor='lightblue']
#   Dementia  [label='Dementia\\n(1)', fillcolor='orange']
#   Death     [label='Death\\n(2)', fillcolor='red']
#   
#   # Переходы
#   NoDementia -> Dementia [label='0→1', color='orange', fontcolor='orange', penwidth=2]
#   NoDementia -> Death    [label='0→2', color='red', fontcolor='red', penwidth=2]
#   Dementia -> Death      [label='1→2', color='red', fontcolor='red', penwidth=2]
#   
#   # Подпись снизу
#   subgraph legend {
#     label='Note: Dementia and Death are absorbing states'
#     fontsize=10
#     style=invis
#   }
# }
# ")

# ggsave(
#   filename = "figures_annex/state_space.pdf",
#   dpi = 300
# )
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

g <- grViz("
digraph dementia_states {

graph [
  rankdir = LR,
  bgcolor = white
]

node [
  shape = circle,
  style = filled,
  fontname = Arial,
  fontsize = 10,
  width = 1.1,
  height = 1.1
]

edge [
  fontname = Arial,
  fontsize = 10,
  penwidth = 1.5
]

Healthy [
  label = 'Dementia-free\\n(0)',
  fillcolor = '#D9EAF7'
]

Dementia [
  label = 'With Dementia\\n(1)',
  fillcolor = '#FCE4D6'
]

Death [
  label = 'Death\\n(2)',
  fillcolor = '#E6E6E6'
]

Healthy -> Dementia [
  label = '0 → 1'
]

Healthy -> Death [
  label = '0 → 2'
]

Dementia -> Death [
  label = '1 → 2'
]

subgraph cluster_legend {
  label = 'Note: Dementia and death are absorbing states.'
  fontname = Arial
  fontsize = 8
  style = invis
}

}
")

# сохранить как PDF
svg <- export_svg(g)

rsvg_pdf(
  charToRaw(svg),
  file = "figures_annex/state_space_diagram1.pdf"
)
# ------------------------------------------------------------------------ #
# ADJ and UNADJ trans hazards
haz_unadj <- read_csv("Data/model2/unadj_haz_replicates.csv.gz", show_col_types = FALSE) %>% 
  filter(year %in% c(2004, 2019)) %>% 
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
  filter(year %in% c(2004, 2019)) %>% 
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

# before and after adjustment hazards all 3
# FOr incidence adjusted and unadjusted are equal
# REMOVE
# FOund another question dementia free hazards are higher for men than dementia in older ages
# haz %>% 
#   full_join(haz_unadj) %>% 
#   mutate(female = ifelse(female == 1, "Men", "Women"),
#          female = factor(female, levels = c("Men", "Women"))) %>% 
#   separate(Transtion, c("From", "To")) %>%
#   mutate(across(c(From, To), ~ case_when(
#     . == "U" ~ "Dementia",
#     . == "H" ~ "Dementia-free",
#     . == "D" ~ "Death",
#   ))) %>% 
#   filter(year %in% c(2019)) %>%
#   filter(To != "Dementia") %>%
#   # filter(To == "Dementia") %>% 
#   # select(-lower, -upper) %>% 
#   # pivot_wider(names_from = Hazards,
#   #             values_from = p_med ) %>% 
#   # select(-year)
#   ggplot(aes(x = age, y = p_med, color = To, linetype = Hazards)) +
#   geom_line(linewidth = 1) + 
#   geom_ribbon(
#     aes(ymin = lower, ymax = upper, fill = To),
#     alpha = 0.12,
#     color = NA
#   ) +
#   scale_fill_brewer(
#     palette = "Dark2"
#   )+
#   scale_color_brewer(
#     palette = "Dark2"
#   ) +
#   theme_bw(base_family = "Arial", base_size = 10) +
#   scale_y_log10() +
#   facet_grid(female ~ From, switch = "y") + 
#   theme(
#     strip.placement = "outside",
#     strip.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.spacing = unit(1.1, "lines"),
#     legend.position = "bottom",
#     axis.title = element_text(size = 10, color = "black", face = "bold"),
#     axis.text = element_text(size = 10, color = "black"),
#     strip.text = element_text(size = 10, face = "bold", color = "black"),
#     legend.title = element_text(size = 10, face = "bold", color = "black"),
#     legend.text = element_text(size = 10, color = "black"),
#     legend.box = "horizontal")+
#   labs(
#     x = "Age (years)",
#     y = "Mortality hazard (log scale), 95% CI",
#     color = "To",
#     fill = "To"
#   )

haz %>% 
  full_join(haz_unadj) %>% 
  mutate(female = ifelse(female == 1, "Men", "Women"),
         female = factor(female, levels = c("Men", "Women"))) %>% 
  separate(Transtion, c("From", "To")) %>%
  mutate(across(c(From, To), ~ case_when(
    . == "U" ~ "Dementia",
    . == "H" ~ "Dementia-free",
    . == "D" ~ "Death",
  ))) %>% 
  filter(year %in% c(2004, 2019)) %>%
  filter(To == "Death") %>%
  select(-To) %>% 
  mutate(year = as.factor(year)) %>% 
  ggplot(aes(x = age, y = p_med, color = year, linetype = Hazards)) +
  geom_line(linewidth = 1) + 
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = year),
    alpha = 0.12,
    color = NA
  ) +
  scale_fill_brewer(
    palette = "Dark2"
  )+
  scale_color_brewer(
    palette = "Dark2"
  ) +
  theme_bw(base_family = "Arial", base_size = 10) +
  scale_y_log10() +
  facet_grid(female ~ From, switch = "y") + 
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.1, "lines"),
    legend.position = "bottom",
    axis.title = element_text(size = 10, color = "black", face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 10, face = "bold", color = "black"),
    legend.title = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.box = "horizontal")+
  labs(
    x = "Age (years)",
    y = "Mortality hazard (log scale), 95% CI",
    color = "Year",
    fill = "Year"
  )

ggsave(
  filename = "figures_annex/haz_2004_2019_adj_uadj.pdf",
  device = cairo_pdf,
  dpi = 300
)

# haz %>% 
#   full_join(haz_unadj) %>% 
#   mutate(female = ifelse(female == 1, "Men", "Women"),
#          female = factor(female, levels = c("Men", "Women"))) %>% 
#   separate(Transtion, c("From", "To")) %>%
#   mutate(across(c(From, To), ~ case_when(
#     . == "U" ~ "Dementia",
#     . == "H" ~ "Dementia-free",
#     . == "D" ~ "Death",
#   ))) %>%
#   filter(To != "Dementia") %>% 
#   mutate(Year = as.factor(year)) %>% 
#   ggplot(aes(x = age, y = p_med, color = Year, linetype = Hazards)) +
#   geom_line(linewidth = 1) + 
#   geom_ribbon(
#     aes(ymin = lower, ymax = upper, fill = Year),
#     alpha = 0.12,
#     color = NA
#   ) +
#   scale_fill_brewer(
#     palette = "Dark2"
#   )+
#   scale_color_brewer(
#     palette = "Dark2"
#   ) +
#   theme_bw(base_size = 14) +
#   scale_y_log10() +
#   facet_grid(female ~ From, switch = "y") + 
#   theme(
#     strip.placement = "outside",
#     strip.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.spacing = unit(1.1, "lines"),
#     axis.title = element_text(size = 12, color = "black"),
#     axis.text = element_text(size = 10, color = "black"),
#     strip.text = element_text(
#       face = "bold",
#       size = 10, 
#       color = "black"
#     ),
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold", color = "black"),
#     legend.box = "horizontal",
#     plot.title = element_text(
#       face = "bold",
#       size = 14, 
#       color = "black"
#     ))+
#   labs(
#     x = "Age",
#     y = "Transition hazard (log scale)",
#     color = "Year",
#     fill = "Year"
#   )
# 
# ggsave(
#   filename = "figures_annex/hazards_years.pdf",
#   dpi = 300
# )

# ------------------------------------------------------------------------ #
# adjusted transtiion probabilities if needed can improve quality
# I think it is not needed
probs <- read_csv("Data/model2/probs.csv.gz")
probs %>% 
  filter(year == 2019) %>% 
  filter(to > from) %>%
  group_by(female, year, age, from, to) %>% 
  summarize(
    p_med = median(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    female = ifelse(female == 1, "Men", "Women"),
    female = factor(female, levels = c("Men", "Women")),
    From = case_when(
      from == 1 ~ "Dementia-free",
      from == 2 ~ "Dementia"
    ),
    To = case_when(
      to == 1 ~ "Dementia-free",
      to == 2 ~ "Dementia",
      to == 3 ~ "Death"
    )
  ) %>%
  ggplot(
    aes(
      x = age,
      y = p_med,
      color = To
    )
  ) +
  geom_line(
    linewidth = 1
  ) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper,
      fill = To
    ),
    alpha = 0.12,
    color = NA
  ) +
  scale_fill_brewer(
    palette = "Dark2"
  ) +
  scale_color_brewer(
    palette = "Dark2"
  ) +
  theme_bw(base_family = "Arial", base_size = 10) +
  scale_y_log10() +
  facet_grid(
    female ~ From,
    switch = "y"
  ) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.1, "lines"),
    axis.title = element_text(size = 10, color = "black", face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 10, face = "bold", color = "black"),
    legend.title = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  labs(
    x = "Age (years)",
    y = "Transition probability (log scale), 95% CI",
    color = "To",
    fill = "To"
  )

ggsave(
  filename = "figures_annex/trans_prob.pdf",
  device = cairo_pdf,
  dpi = 300
)



prop <- read_csv("Data/model2/prop50.csv.gz")

# prop |>
#   mutate(
#     female = factor(
#       female,
#       levels = c(0, 1),
#       labels = c("Men", "Women")
#     )
#   ) |>
#   mutate(across(c(prop50, lower, upper), ~ 1 - .)) %>% 
#   filter(year %in% c(2004, 2019))


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
    y = "Proportion of remaining life spent with dementia at age 50, \n with 95% CI",
    color = "Gender",
    fill = "Gender"
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  theme_bw(base_family = "Arial", base_size = 10) +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10, color = "black", face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 10, face = "bold", color = "black"),
    legend.title = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    panel.spacing = unit(1.2, "lines")
  )


ggsave(
  filename = "figures_annex/proportion.pdf",
  device = cairo_pdf,
  dpi = 300
)

# ------------------------------------------------------------------- #
# Prevalence weighted mortality before and after ADJ.
# CHECK!
prev <- read_csv("Data/model2/prev_replicates.csv.gz",show_col_types = FALSE)

hazards_unadj <- read_csv("Data/model2/unadj_haz_replicates.csv.gz", show_col_types = FALSE) %>% 
  filter(to > from) |> # create c from and to columns
  mutate(from = ifelse(from == 1, "H", "U"),
         to   = case_when(
           to == 1 ~ "H",
           to == 2 ~ "U",
           to == 3 ~ "D"
         )) %>% 
  mutate(transition = paste0(from, "_", to)) %>%
  select(-from, -to) %>%
  pivot_wider(names_from = transition, values_from = haz) %>%
  left_join(
    prev,
    by = join_by(replicate, age, female, year)
  ) %>% 
  select(-H_U) %>%
  mutate(
    mx = (1 - prevalence) * H_D + prevalence  * U_D) %>%
  arrange(replicate, female, year, age) 

# Absolutely the same procedure just for adjusted
hazards_adj <- read_csv("Data/model2/adj_haz_replicates.csv.gz", show_col_types = FALSE) %>%
  filter(to > from) |> # create c from and to columns
  mutate(from = ifelse(from == 1, "H", "U"),
         to   = case_when(
           to == 1 ~ "H",
           to == 2 ~ "U",
           to == 3 ~ "D"
         )) %>% 
  mutate(transition = paste0(from, "_", to)) %>%
  select(-from, -to) %>%
  pivot_wider(names_from = transition, values_from = haz) %>%
  left_join(
    prev,
    by = join_by(replicate, age, female, year)
  ) %>% 
  select(-H_U) %>%
  mutate(
    mx = (1 - prevalence) * H_D + prevalence  * U_D) %>%
  arrange(replicate, female, year, age) 


unique(diff(sort(unique(hazards_unadj$age))))

summary(hazards_unadj$mx)

hazards_unadj %>%
  group_by(replicate, year, female) %>%
  summarise(n=n())

hazards_unadj %>%
  group_by(year,female) %>%
  summarise(first=min(age), last=max(age))

# unite hazards and calculate overall mx
mx_plot <- bind_rows(
  hazards_unadj %>% 
    group_by(female, year, age) %>%
    summarise(
      p_med = median(mx),
      lower = quantile(mx, 0.025),
      upper = quantile(mx, 0.975),
      .groups = "drop"
    ) %>% 
    mutate(Hazards = "Before adjustment"),
  hazards_adj   %>%
    # summarise bands     
    group_by(female, year, age) %>%
    summarise(
      p_med = median(mx),
      lower = quantile(mx, 0.025),
      upper = quantile(mx, 0.975),
      .groups = "drop"
    ) %>% 
    mutate(Hazards = "After adjustment")
)
# ------------------------------------------------------------------- #
# Resulting e50 calculation before and after
unadj_e50 <- hazards_unadj %>% 
  group_by(replicate, year, female) %>%
  summarise(
    # I assume quarters     
    e50 = sum(exp(-cumsum(mx * 0.25))) * 0.25,
    .groups = "drop"
  )

unadj_e50_summary <- unadj_e50 %>%
  group_by(female, year) %>%
  summarise(
    e50 = median(e50, na.rm = TRUE),
    lower = quantile(e50, 0.025, na.rm = TRUE),
    upper = quantile(e50, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# absolutely the same for adjusted
adj_e50 <- hazards_adj %>% 
  group_by(replicate, year, female) %>%
  summarise(
    e50 = sum(exp(-cumsum(mx * 0.25))) * 0.25,
    .groups = "drop"
  )

adj_e50_summary <- adj_e50 %>%
  group_by(female, year) %>%
  summarise(
    e50 = median(e50, na.rm = TRUE),
    lower = quantile(e50, 0.025, na.rm = TRUE),
    upper = quantile(e50, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# join
e50_plot <- bind_rows(
  unadj_e50_summary %>% 
    mutate(Hazards = "Before adjustment"),
  adj_e50_summary %>%
    mutate(Hazards = "After adjustment")
) %>%
  mutate(female = ifelse(female == 1, "Female", "Male"))

# Something wrong
e50_plot %>% 
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
