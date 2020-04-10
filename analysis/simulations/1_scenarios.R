##******************************************##
##     Scenario building for simulations    ##
##******************************************##



# Parameters varied -------------------------------------------------------

# Number of replications per scenario, based on monte-carlo error bias
n_sim <- 160 


# All possible parameters to vary
n <- c(n = 500, "large" = 2000)
prop_miss <- c("low" = .10, "high" = .50)
beta1 <- c("null" = 0, "med" = 0.5, "large" = 1)
mech <- c("MCAR", "MAR", "MNAR", "MAR_GEN")
X_level <- c("continous", "binary")
rho <- c("weak" = 0.1, "strong" = 0.5)
eta1 <- c("weaker" = -1, "strong" = -2) # To incorporate MAR_GEN
haz_shape <- c("similar", "different")



# Re-visit 07/04/2019 -----------------------------------------------------


# Make the full factorial and label the pilots
full_factorial <- expand.grid("n" = n,
                              "prop_miss" = prop_miss,
                              "beta1" = beta1,
                              "miss_mech" = mech,
                              "X_level" = X_level,
                              "rho" = rho,
                              "eta1" = eta1,
                              "haz_shape" = haz_shape) %>% 
                         
  # Eta 1 not relevant for MCAR. removes duplicates
  dplyr::mutate(eta1 = ifelse(miss_mech == "MCAR", NA, eta1)) %>% 
  dplyr::distinct() %>% 

  # Pick out the pilot scenarios
  dplyr::mutate(
    pilot = dplyr::case_when(
      beta1 == 0.5 ~ 0,
      n == 500 ~ 0,
      rho == 0.1 ~ 0,
      X_level == "binary" ~ 0, 
      prop_miss == 0.1 ~ 0,
      haz_shape == "different" ~ 0,
      TRUE ~ 1
    )
  )  
  
# Subset pilots and given them a scen number
pilot_scenarios <- full_factorial %>% 
  
  dplyr::filter(pilot == 1) %>% 
  
  # We run the pilot ones first before deciding which other scenarios to run
  dplyr::arrange(dplyr::desc(pilot)) %>% 
  dplyr::filter(pilot == 1) %>% 
  dplyr::mutate(scen_num = 1:dplyr::n())


# Decisions after running these:
# - beta1 = 0 not of primary interest
# - m = 50 is enough for imputations
# - rho = 0.5 kept

# After running two scenarios with n = 500
# - Only run n = 500 for pilot scens, as check 
# - use more MC reps, for interpretation (needs to be check)


# Set up copy of pilots with n = 500
# No scen number yet as we do not yet no n_reps
pilots_n500 <- pilot_scenarios %>% 
  dplyr::mutate(n = 500) %>% 
  dplyr::select(-scen_num)


# Now set-up the rest 
scenarios_remaining <- full_factorial %>% 
  
  # Remove pilots
  dplyr::filter(pilot == 0)  %>% 
  
  # Remove other criteria
  dplyr::filter(
    rho == 0.5,
    beta1 != 0,
    n == 2000
  ) %>% 
  
  # Sort by X_level and haz_shape (for batching)
  dplyr::arrange(X_level, haz_shape) %>% 

  # Add scen_num
  dplyr::mutate(scen_num = 1:dplyr::n() + max(pilot_scenarios$scen_num))


# Bind it all together (n = 500 will come later)
scenarios_updated <- dplyr::bind_rows(
  pilot_scenarios,
  scenarios_remaining
) %>% 
  dplyr::mutate(seed = scen_num * n_sim)


# Save as .RDS file to load in
saveRDS(scenarios_updated, file = "inst/testdata/scenarios.rds")

# Need n = 500 scens last...


# Other -------------------------------------------------------------------



# This is for using summarise sims on subset
# Replace 1 and 2 by some actual number
list.files(
  path = "./analysis/sim_results/summarised_reps/estimates/",
  pattern = ".*(scen28_)|(scen57_).*.rds", 
  full.names = T
)