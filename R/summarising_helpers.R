##******************************************##
## Functions for summarising after sims run ##
##******************************************##



# Test function for rmse and emp se
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



#bind_rows(final_test) %>%
#  group_by(var, analy, m) %>%
# mutate(nsim = n()) %>% 
# mutate(rmse = rmse(coef_i1, true, nsim),
#        emp_se = emp_SE(coef_i1, mean(coef), nsim)) %>% 
# mutate(se_emp = sd(coef)) %>%
#  summarise_all(~ round(mean(.), 3)) 

# After running 
# print(xtable(manip), digits = 3), include.rownames = FALSE)



# Loading RDS files -------------------------------------------------------


# Get names of files
#list_files <- list.files(".", pattern = 'results', full.names = T)


# Function coerces them into a list
loadRDS <- function(fileName){
  #loads an RDS file, and returns it
  readRDS(fileName)
  get(ls()[ls() != "fileName"])
}


