##******************************************##
##     Scenario building for simulations    ##
##******************************************##



# Parameters varied -------------------------------------------------------


# Number of replications per scenario, based on monte-carlo error bias
n_sim <- 160 

# Parameters to vary
n <- c("small" = 500, "large" = 2000)
prop_miss <- c("low" = .10, "high" = .50)
beta <- c("null" = 0, "med" = 0.5, "large" = 1)
mech <- c("MCAR", "MAR", "MNAR", "MAR_GEN")
X_level <- c("continous", "binary")
rho <- c("weak" = 0.1, "strong" = 0.5)
eta1 <- c("little" = -1, "strong" = -2) # To incorporate MAR_GEN


# Make a scenarios grid, this is full factorial case
scenarios <- expand.grid("n" = n,
                         "prop_miss" = prop_miss,
                         "beta1" = beta,
                         "miss_mech" = mech,
                         "X_level" = X_level,
                         "rho" = rho,
                         "eta1" = eta1) %>% 
  
  # Eta 1 not relevant for MCAR. removes duplicates
  dplyr::mutate(eta1 = ifelse(miss_mech == "MCAR", NA, eta1)) %>% 
  dplyr::distinct() %>% 
  
  # The pilot 14 scenarios, by *excluding* the following levels
  dplyr::mutate(
    pilot = ifelse(
      n == 500 | beta == 0.5 | X_level == "binary" |
        prop_miss == 0.1 | rho == 0.1,
      0, 1
    )
  ) %>% 
  
  # We run the pilot ones first before deciding which other scenarios to run
  dplyr::arrange(dplyr::desc(pilot)) %>% 
  
  # Assign scenario number and base seed
  dplyr::mutate(scen_num = 1:dplyr::n(),
                seed = n_sim * scen_num)

rm(n, prop_miss, beta, mech, X_level, rho, eta1)
  
# This is full factorial
scenarios

# 
pilot_scens <- scenarios[scenarios$pilot == 1, ]

# Now take subset after discussion 25/03/2020
scenarios_sub <- scenarios %>% 
  dplyr::filter(
    !pilot,
    rho == 0.5,
    beta1 != 0
  ) %>% 
  dplyr::bind_rows(pilot_scens) %>% 
  dplyr::arrange(dplyr::desc(pilot)) %>% 
  
  # Re-do the indices 
  dplyr::mutate(
    scen_num = 1:dplyr::n(),
    seed = n_sim * scen_num
  )


# Save as .RDS file to load in
saveRDS(scenarios_sub, file = "inst/testdata/scenarios.rds")


# This is for using summarise sims on subset
# Replace 1 and 2 by some actual number
list.files(
       path = "./analysis/sim_results/summarised_reps/estimates/",
       pattern = ".*(scen28_)|(scen57_).*.rds", 
       full.names = T
   )
