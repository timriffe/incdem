# In this file I will try to check if the patern of transition hazards 
# for H-D and U-D holds in other data.


# In their results, mortality decreases for both men and women, 
# and for people both with and without dementia. 
#
# In our case, mortality among healthy individuals decreases for men 
# but increases slightly for women, especially at older ages. 
# Mortality among people with dementia increases for men and 
# decreases slightly for women, particularly at younger ages. 
# This difference is consistent across Models 2 and 3.



# I start by analysing the HMD all cause USA mortality.
# HMD should mostly map to HD. 

# Compare 2004 mortality rates with 2019
# HMD overall mortality rate declines for both sexes
# in our case HD declines for M and increase in F
mx_usa <- read_table("Data/hmd/Mx_1x1.txt", skip = 1) %>% 
  select(-Total) %>% 
  filter(Year %in% c(2004:2019))

mx_usa %>%
  pivot_longer(c(Female, Male),
               names_to  = "sex",
               values_to = "mx") %>% 
  filter(Year %in% c(2004, 2019), Age > 49) %>%
  mutate(Age = as.numeric(Age)) %>%
  mutate(Year = as.factor(Year)) %>% 
  ggplot(aes(x = Age, y = mx, group = Year, color = Year)) + 
  geom_line() + 
  scale_y_log10() +
  facet_wrap(~ sex)  + 
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave("hmd_log_hazards.jpeg")


# Lets now check the mortality from the human cause of death database
# They have a full lost of causes that has dementia
# they recently updated the ddata so we are golden
# Vascular dementia F01 and Unspecified dementia F03
mx_usa_c <- read_csv("Data/hmd/USA_m_full_idr (1).csv") %>%
  filter(year %in% c(2004, 2019)) %>% 
  select(-c(total, country, list, agf)) %>% 
  pivot_longer(-c(year, sex, cause),
               names_to  = "age",
               values_to = "mx") %>% 
  mutate(age = parse_number(age)) %>% 
  filter(age > 49) %>% 
  filter(cause %in% c("F01", "F03")) %>% 
  filter(sex != 3) %>% 
  mutate(sex = ifelse(sex == 1, "male", "female"))

# Mortality in people with dementia increase for both dementias
# in our case decline in F and increase in M
mx_usa_c %>% 
  mutate(year = as.factor(year)) %>%
  ggplot(aes(x = age, y = mx, group = year, color = year)) +
  geom_line() +
  scale_y_log10() +
  facet_grid(sex ~ cause) +
  theme_bw() +
  theme(legend.position = "bottom")



# library("LEdecomp")
# 
# decomp_usa <- LEdecomp::US_data_CoD
# 
# decomp_usa %>% 
#   count(cause_id)
# 
# ggsave("trans_hazards.jpeg")
# ggsave("log_mort_rate_usa_hmd.jpeg")




