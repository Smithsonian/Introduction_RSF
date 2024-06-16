setwd("C:/Users/stabachj/Dropbox (Smithsonian)/TheContinental_VARIOGRAMS")

library(stringr)
library(tidyverse)

# Search the list
all <- list.files(path='.', pattern='png', all.files=FALSE, full.names=FALSE)
all <- as_tibble(as.data.frame(all))

all <-
  all %>% 
  rename(name = all) %>% 
  mutate(
    animal = str_sub(name,
            start = 12, 
            end = -5),
    resident = 0,
    completed = 0) %>% 
  write_csv(file = "Resident_List.csv")
