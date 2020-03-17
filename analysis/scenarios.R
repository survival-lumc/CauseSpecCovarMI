##******************************************##
##     Scenario building for simulations    ##
##******************************************##



# Parameters varied -------------------------------------------------------


n_sim <- 160 # change this based on a monte carlo error!!!



#n <- c("small" = 500, "large" = 2000)
n <- c("large" = 2000)
#prop_miss <- c("low" = .10, "high" = .5)
prop_miss <- c("high" = .5)

#beta <- c("null" = 0, "med" = 0.5, "large" = 1)
beta <- c("null" = 0, "large" = 1)
mech <- c("MCAR", "MAR", "MNAR", "MAR-GEN")
#X_level <- c("continous", "binary")
X_level <- c("continous")
#rho <- c("weak" = 0.1, "strong" = 0.5)
rho <- c("strong" = 0.5)
eta1 <- c("little" = 1, "strong" = 2)


# The numbers
scens <- expand.grid(
  "n" = n,
  "prop_miss" = prop_miss,
  "beta1" = beta,
  "miss_mech" = mech,
  "X_level" = X_level,
  "rho" = rho,
  "eta1" = eta1
) %>% 
  mutate(eta1 = ifelse(miss_mech == "MCAR", NA, eta1),
         scen_num = 1:dplyr::n(),
         seed = n_sim * scen_num) %>% 
  dplyr::distinct()

rm(n_sim, n, prop_miss, beta, mech, X_level, rho, eta1)
  
# With the labels
scens %>% 
  map_if(is.numeric, names) %>% 
  as_tibble() 


# Save as .RDS file to load in
# saveRDS(scenarios, file = "scenarios.rds")