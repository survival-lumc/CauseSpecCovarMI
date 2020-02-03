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

# Add functions empse and rmse
emp_SE <- function(theta_i1, theta_hat, nsim) {
  
  sq_diffs <- (theta_i1 - theta_hat)^2 
  term <- sum(sq_diffs) / (nsim - 1)
  return(sqrt(term))
}

rmse <- function(theta_i1, true, nsim) {
  
  sq_diffs <- (theta_i1 - true)^2 
  term <- sum(sq_diffs) / (nsim - 1)
  return(sqrt(term))
}


# Bind list together
final_dat <- bind_rows(lapply(list_files, loadRData)) %>% 
  group_by(scen_name, var, analy, m) %>% 
  mutate(nsim = n()) %>% 
  mutate(rmse = rmse(coef_i1, true, nsim),
         emp_se = emp_SE(coef_i1, mean(coef), nsim)) %>% 
  mutate(se_emp = sd(coef)) %>%
  summarise_all(~ round(mean(.), 3)) %>%
  as.data.frame()


# Save in chosen directory
save(final_dat, file = "/exports/molepi/users/efbonneville/mi_comprisks/final_dat.RData")