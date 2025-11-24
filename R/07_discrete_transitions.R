
source("R/00_dependencies.R")
source("R/01_functions.R")



# ----------------------------------------------
# model 1, 
# ----------------------------------------------
haz <- read_csv("Data/model1/adj_haz_replicates.csv.gz") |> 
  mutate(from = substr(transition,2,2) |> as.integer(),
         to = substr(transition,3,3) |> as.integer()) |> 
  select(-transition)
object.size(haz) |> print(units="Mb")

probs1f <- 
  haz |> 
  filter(female == 1) |> 
  hazards_to_discrete(
    age_interval = 0.25,
    id_cols      = c("age", "replicate", "female", "period5"),
    n_cores      = 7,
    parallel     = "mclapply")

write_csv(probs1f, file = "Data/model1/probs_f.csv.gz")
rm(probs1f);gc()

probs1m <- 
  haz |> 
  filter(female == 0) |> 
  hazards_to_discrete(
    age_interval = 0.25,
    id_cols      = c("age", "replicate", "female", "period5"),
    n_cores      = 6,
    parallel     = "mclapply")

write_csv(probs1m, file = "Data/model1/probs_m.csv.gz")
rm(probs1m);gc()

probs <- vroom(c("Data/model1/probs_m.csv.gz","Data/model1/probs_f.csv.gz"))
write_csv(probs, "Data/model1/probs.csv.gz")
rm(probs);gc()
unlink("Data/model1/probs_m.csv.gz");unlink("Data/model1/probs_f.csv.gz")
# ----------------------------------------------
# model 2, loop over years and sex:
# ----------------------------------------------
haz <- read_csv("Data/model2/adj_haz_replicates.csv.gz")
object.size(haz) |> print(units="Mb")

yrs <- haz$year |> unique() |> sort()
sxs <- c(0,1)
for (y in yrs){
  for (s in sxs){
  hazy <- haz |> filter(year == y & female == s)
  haz <- haz |> filter(!(year == y & female == s));gc()
  if (nrow(hazy)>0){
  probsy <- hazards_to_discrete(
    hazy,
    age_interval = 0.25,
    id_cols      = c("replicate","female","year","age"),
    n_cores      = 4,
    parallel     = "mclapply")
  rm(hazy);gc()
  namey <- paste0("probs_",y,"_",s,".csv.gz")
  write_csv(probsy, file = file.path("Data/model2/tmp",namey))
  rm(probsy);gc()
  }
}}
rm(haz);gc()
probs_files<- file.path("Data/model2/tmp",dir("Data/model2/tmp"))
probs <- vroom(probs_files)
write_csv(probs,file = "Data/model2/probs.csv.gz")
rm(probs);gc()
unlink("Data/model2/tmp", recursive=TRUE)

# ----------------------------------------------
# model 3, loop over years and sex:
# ----------------------------------------------

haz <- read_csv("Data/model3/adj_haz_replicates.csv.gz")
object.size(haz) |> print(units="Mb")
dir.create("Data/model3/tmp")
yrs <- haz$year |> unique() |> sort()
sxs <- c(0,1)
for (y in yrs){
  for (s in sxs){
    hazy <- haz |> filter(year == y & female == s)
    haz <- haz |> filter(!(year == y & female == s));gc()
    if (nrow(hazy)>0){
      probsy <- hazards_to_discrete(
        hazy,
        age_interval = 0.25,
        id_cols      = c("replicate","female","year","age"),
        n_cores      = 4,
        parallel     = "mclapply")
      rm(hazy);gc()
      namey <- paste0("probs_",y,"_",s,".csv.gz")
      write_csv(probsy, file = file.path("Data/model3/tmp",namey))
      rm(probsy);gc()
    }
  }}
rm(haz);gc()
probs_files<- file.path("Data/model3/tmp",dir("Data/model3/tmp"))
probs <- vroom(probs_files)
write_csv(probs,file = "Data/model3/probs.csv.gz")
rm(probs);gc()
unlink("Data/model3/tmp", recursive=TRUE)


