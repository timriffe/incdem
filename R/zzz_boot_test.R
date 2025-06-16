library(rsample)
library(tidvyerse)
# maybe pre-define which ids enter each boot replicate, using a minimal
# data structure

set.seed(123)
df <- data.frame(
  id = rep(1:100, each = 4),
  wave = rep(1:4, 100),
  weight = runif(100)[rep(1:100, each = 4)],
  y = rnorm(400)
)

source("R/zzz_boot_test_functions.R")
boot2 <- group_bootstraps2(df, group = "id", weight = "weight", times = 100)

# View some results
id_freq <- purrr::map_dfr(boot2$splits, function(split) {
  tibble(id = analysis(split)$id)
}) %>%
  count(id) %>%
  arrange(desc(n))

id_freq <- left_join(id_freq, df %>% distinct(id, weight), by = "id")|> 
  mutate(np = n / sum(n),
         weightp = weight / sum(weight))
id_freq |> 
  ggplot(aes(x=np,y=weightp)) +
  coord_equal() +
  geom_point() +
  geom_abline(slope=1,intercept=0)


