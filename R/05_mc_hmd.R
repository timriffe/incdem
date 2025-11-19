
library(HMDHFDplus)
library(janitor)

Deaths <- readHMDweb("USA","Deaths_1x1", 
                 username = Sys.getenv("pw"), 
                 password = Sys.getenv("pw")) |> 
  clean_names() |> 
  select(-total) |> 
  filter(between(age, 40, 110),   # these could be global variables set elsewhere
         between(year, 2004, 2019)) |> 
  pivot_longer(c(male, female),names_to = "sex", values_to = "Dx")

Exposures <- readHMDweb("USA","Exposures_1x1", 
                 username = Sys.getenv("pw"), 
                 password = Sys.getenv("pw"))|> 
  clean_names() |> 
  select(-total) |> 
  filter(between(age, 40, 110),   # these could be global variables set elsewhere
         between(year, 2004, 2019)) |> 
  pivot_longer(c(male, female), names_to = "sex", values_to = "exposure")


# make a helper that smooths mortality ages 40-110,
# interpolating to .25 age intervals (.25 should also be a global parameter)

# when done, select down to ages 50-100










