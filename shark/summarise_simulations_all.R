#!/bin/bash


#####################################
## Reads in all .Rdata simulations ##
##       and summarises them       ##
#####################################


# Saves everything
library(tidyverse)


# Get names of files
list_files <- list.files(".", pattern = 'results', full.names = T)


# Function coerces them into a list
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# Bind list together
final_dat <- bind_rows(lapply(list_files, loadRData)) %>% 
  group_by(scen_name, var, analy, m) %>% 
  mutate(se_emp = sd(coef)) %>% 
  summarise_all(~ round(mean(.), 3)) %>% 
  ungroup() %>% 
  as.data.frame()


# Save in chosen directory
save(final_dat, file = "/exports/molepi/users/efbonneville/mi_comprisks/final_dat.RData")