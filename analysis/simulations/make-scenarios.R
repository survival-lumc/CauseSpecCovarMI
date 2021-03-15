##******************************************##
##     Scenario building for simulations    ##
##******************************************##


library(CauseSpecCovarMI)


# Parameters varied -------------------------------------------------------


# Number of replications per scenario, based on monte-carlo error bias
n_sim <- 160 

# All possible parameters to vary
n <- c(n = 500, "large" = 2000)
prop_miss <- c("low" = .10, "high" = .50)
beta1 <- c("null" = 0, "med" = 0.5, "large" = 1)
mech <- c("MCAR", "MAR", "MNAR", "MAR_GEN")
X_level <- c("continuous", "binary")
rho <- c("weak" = 0.1, "strong" = 0.5)
eta1 <- c("weaker" = -1, "strong" = -2) # To incorporate MAR-T
haz_shape <- c("similar", "different")


# Creating full factorial, and pilot scenarios ----------------------------


# Make the full factorial and label the pilots
full_factorial <- expand.grid(
  "n" = n,
  "prop_miss" = prop_miss,
  "beta1" = beta1,
  "miss_mech" = mech,
  "X_level" = X_level,
  "rho" = rho,
  "eta1" = eta1,
  "haz_shape" = haz_shape
) %>% 
  data.table() %>% 
  
  # Eta1 not relevant for MCAR, remove resulting duplicates
  .[, eta1 := ifelse(miss_mech == "MCAR", NA, eta1)] %>% 
  unique() %>% 

  # Pick out the pilot scenarios
  .[, pilot := data.table::fcase(
    beta1 == 0.5, 0,
    n == 500, 0,
    rho == 0.1, 0,
    X_level == "binary", 0, 
    prop_miss == 0.1, 0,
    haz_shape == "different", 0,
    default = 1
  )]

# Subset pilots and given them a scenario number
pilot_scenarios <- data.table::copy(full_factorial) %>% 
  .[pilot == 1] %>% 
  .[order(pilot), scen_num := 1:.N]

# Decisions after running these:
# - beta1 = 0 not of primary interest
# - m = 50 is enough for imputations
# - rho = 0.5 kept

# After running two scenarios with n = 500
# - Only run n = 500 for pilot scenarios, as check 


# All scenarios -----------------------------------------------------------


# Remaining scenarios
scenarios_remaining <- data.table::copy(full_factorial) %>% 
  
  # Remove pilots and other criteria
  .[pilot == 0 & rho == 0.5 & beta1 != 0 & n == 2000] %>% 
  .[order(X_level, prop_miss, haz_shape), scen_num := 1:.N + max(
    pilot_scenarios$scen_num
  )]

# Set up copy of pilots with n = 500
n_sim_500 <- 400 # 0.2^2/0.01^2
pilots_n500 <- data.table::copy(pilot_scenarios) %>% 
  .[, ':=' (
    n = 500,
    scen_num = 1:.N + max(scenarios_remaining$scen_num)
  )] 
  
# Bind it all together 
scenarios_full <- rbind(
  pilot_scenarios, 
  scenarios_remaining,
  pilots_n500
) %>% 
  .[n == 2000, seed := scen_num * n_sim] %>% 
  .[n == 500, seed := scen_num * n_sim_500] %>% 
  .[order(scen_num)]

# Note we do not report results of n = 500

# Save as .RDS file to load in
saveRDS(scenarios_full, file = "data/scenarios.rds")
