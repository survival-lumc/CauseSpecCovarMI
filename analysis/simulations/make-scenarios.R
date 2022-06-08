##******************************************##
##     Scenario building for simulations    ##
##******************************************##


# This file contains how the scenarios grid was originally prepared for the 
# simulation study. It contains the 112 scenarios described in the article,
# plus an additional 7 scenarios with beta_1 = 0 (to check the null).


# Parameters varied/fixed -------------------------------------------------


n_sim <- 160 # Number of replications per scenario, based on monte carlo error bias
n <- c("large" = 2000) # Sample size 
prop_miss <- c("low" = .10, "high" = .50) # Proportion missing data
beta1 <- c("null" = 0, "med" = 0.5, "large" = 1) # Coefficient on X in model for REL
mech <- c("MCAR", "MAR", "MNAR", "MAR_GEN") # Missingness mechanism
X_level <- c("continuous", "binary") # X covariate type
rho <- c("strong" = 0.5) # Correlation between X and Z
eta1 <- c("weaker" = -1, "strong" = -2) # "Strength" of the missingness mechanism
haz_shape <- c("similar", "different") # Form of the hazard shapes


# Creating full factorial, and pilot scenarios ----------------------------


# Make the full factorial and label the pilots
full_factorial <- data.table::data.table(
  expand.grid(
    "n" = n,
    "prop_miss" = prop_miss,
    "beta1" = beta1,
    "miss_mech" = mech,
    "X_level" = X_level,
    "rho" = rho,
    "eta1" = eta1,
    "haz_shape" = haz_shape
  )
) 
  
# Eta1 not relevant for MCAR, remove resulting duplicates
full_factorial <- unique(full_factorial[, eta1 := ifelse(miss_mech == "MCAR", NA, eta1)])

# Pick out the pilot scenarios
full_factorial[, pilot := data.table::fcase(
  beta1 == 0.5, 0,
  n == 500, 0,
  rho == 0.1, 0,
  X_level == "binary", 0, 
  prop_miss == 0.1, 0,
  haz_shape == "different", 0,
  default = 1
)]

# NOTE: The "pilot" scenarios were named as such since they were used in the 
# testing phase of the simulation study.

# Subset pilots and given them a scenario number
pilot_scenarios <- data.table::copy(full_factorial[pilot == 1])
data.table::setorder(pilot_scenarios, pilot)
pilot_scenarios[, scen_num := 1:.N]

# Decisions after running these:
# - beta1 = 0 not of primary interest
# - m = 50 is enough for imputations
# - rho = 0.5 kept

# After running two scenarios with n = 500
# - Only run n = 500 for pilot scenarios, as check 


# All scenarios -----------------------------------------------------------


# Remaining scenarios
scenarios_remaining <- data.table::copy(full_factorial) 
scenarios_remaining <- scenarios_remaining[pilot == 0 & beta1 != 0]
data.table::setorder(scenarios_remaining, X_level, prop_miss, haz_shape)
scenarios_remaining[, scen_num := 1:.N + max(pilot_scenarios$scen_num)]

# Bind it all together 
scenarios_full <- rbind(pilot_scenarios, scenarios_remaining) 
scenarios_full[, seed := scen_num * n_sim]
data.table::setorder(scenarios_full, scen_num)

# Save as .Rdata to load-in 
#save(scenarios_full, file = "data/scenarios.RData")
