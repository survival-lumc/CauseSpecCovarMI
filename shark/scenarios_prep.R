########################################
## Setting up data frame of scenarios ##
########################################


# library(tidyverse)

# Varied parameters
betas <- c(0, 0.5, 1) # for beta from hazard 1 model
miss <- c("high", "low") # 50 or 10 % missing
mech <- c("MCAR", "MAR_MRR", "MAR_W") # Type of mechanism


# All combinations of factor variables
scenarios <- expand.grid(beta1 = betas, miss = miss, mech = mech) %>% 
  mutate(scen = 1:nrow(.))


# For picking what dataset to simulate for what scenario
sim_scenario <- function(beta1, miss, mech, scen) {
  
  function(seed) {
    
    # Ensure seed is different for each scenario
    set.seed(seed + scen)
    
    # Set missingness
    p <- ifelse(miss == "high", 0.5, 0.1)
    
    # Generate according to mechanism
    if (mech == "MCAR") {
      
      dat <- dat_gener_MRR(N = 500, 
                           X_type = "contin",
                           mus = c(0, 0), 
                           covmat = matrix(c(1, 0.25, 
                                             0.25, 1), nrow = 2), 
                           mech = "MCAR", 
                           p = p, 
                           cause2 = "weib", 
                           vals_t1 = c(0.3, 1, c(beta1, -0.5)), 
                           vals_t2 = c(1.7, 0.5, -0.5, 0.5))
        
    } else if (mech == "MAR_MRR") {
      
      dat <- dat_gener_MRR(N = 500, 
                           X_type = "contin",
                           mus = c(0, 0), 
                           covmat = matrix(c(1, 0.25, 
                                             0.25, 1), nrow = 2), 
                           mech = "MAR", 
                           p = p, 
                           cause2 = "weib", 
                           vals_t1 = c(0.3, 1, c(beta1, -0.5)), 
                           vals_t2 = c(1.7, 0.5, -0.5, 0.5))
      
    } else {
      
      dat <- dat_gener_W(N = 500, 
                         X_type = "contin",
                         mus = c(0, 0), 
                         covmat = matrix(c(1, 0.25, 
                                           0.25, 1), nrow = 2), 
                         mech = "MAR", 
                         pars_MAR = c(1, 1),
                         p = p, 
                         cause2 = "weib", 
                         vals_t1 = c(0.3, 1, c(beta1, -0.5)), 
                         vals_t2 = c(1.7, 0.5, -0.5, 0.5)) 
    }
  }
}


# Add these as functions 
sim_functions <- split(scenarios, 1:nrow(scenarios)) %>% 
  lapply(., function(arguments) sim_scenario(beta1 = arguments$beta1,
                                             miss = arguments$miss,
                                             mech = arguments$mech,
                                             scen = arguments$scen))
  

scenarios <- scenarios %>% 
  mutate(scen_name = paste(as.character(beta1), miss, mech, sep = "_"),
         sim_function = sim_functions)


# Set true values - based on file "truevals_MAR_W.R"
true_betas <- data.frame(t(apply(scenarios, 1, function(row) {
  
  if (row$mech == "MAR_W") {
    if (row$beta1 == 0) {
      
      c(0.017, -0.398)
      
    } else if (row$beta1 == 0.5) {
      
      c(0.393, -0.392)
      
    } else 
      
      c(0.766, -0.388)
    
  } else {
    c(row$beta1, -0.5)
  }
})))


# Append
scenarios <- scenarios %>% 
  mutate(true_b1 = true_betas$X1,
        true_b2 = true_betas$X2)

rm(betas, mech, miss, true_betas)

# Save as .Rdata file to load in
save(scenarios, file = "scenarios.Rdata")


