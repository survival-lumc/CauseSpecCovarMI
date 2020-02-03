#!/bin/bash

##############################
##    Main simulation file  ##
##  this is the one we qsub ##
##############################


# Uncomment and run below if editing file from gitlab directory
# setwd("shark")

# Clear environment
rm(list = ls())

# Read in command line argument
args <- commandArgs(TRUE)
seed <- as.numeric(args[1])


# Load necessary packages - needs loading in this order
pacman::p_load(MASS, tidyverse, mice,
               mitools, survival, smcfcs)

# Load scenarios
load("scenarios.Rdata")

# Data generating functions
source("data_generation/dat_generation_MRR.R")
source("data_generation/dat_generation_W.R")

# Support functions - added true_betas as argument
source("support_functions/mice_pool_diffm.R")
source("support_functions/smcfcs_pool_diffm.R")
source("support_functions/quiet_printcat.R")


# Apply all scenarios for one seed (18 different datasets for 18 scenarios)
res_seed <- lapply(1:nrow(scenarios), function(row) {
  
  # Read in seed locally
  seed <- seed
  
  # Read in parameter of use 
  beta1 <- scenarios$beta1[row]
  mechan <- scenarios$mech[row]
  
  # Generate dataset according to scenario
  dat <- scenarios[row, ]$sim_function[[1]](seed) # this is where the seed comes in
  
  # Set some global parameters
  
  # Reference is different for W
  true_betas <- c(scenarios$true_b1[row],
                  scenarios$true_b2[row])
  
  # Set up MICE matrices
  mpred_ch1 <- mpred_ch12 <- mpred_ch12_int <-  matrix(0, ncol(dat), ncol(dat),
                                                       dimnames = list(names(dat), names(dat)))
  
  ## Ch1 model: MI with Z, eps (as factor) and H1(t) as covars
  mpred_ch1["X1", c("X2", "ev1", "H1")] <- 1
  
  ## Ch12 model: MI with Z, eps, H1(t) & H2(t) as covars
  mpred_ch12["X1", c("X2", "eps", "H1", "H2")] <- 1
  
  # Ch12 model with addition of interactions H x Z
  mpred_ch12_int["X1", c("X2", "eps", "H1", "H2", "H1_X2", "H2_X2")] <- 1
  
  
  #  Global parameters 
  n_rep <- 3 # number of repetitions for mice()
  #m <- c(1, 10, 20) # imputations we are interested in
  m <- c(1,2,3)
  method <- "norm" # continuous case, so use Bayesian linear regression

  
  # Run all methods in each scenario:
  
  
  # Start repetitions for MI methods
  repetitions <- lapply(1:n_rep, function(j) {
    
    # Run imputations, run cause-specific cox, and pool - 
    imp_ch1 <- mice(dat, m = m[length(m)],
                    method = method, 
                    predictorMatrix = mpred_ch1,
                    print = FALSE)
    
    results_ch1 <- mice_pool_diffm(imps = imp_ch1,
                                   m = m, j = j, label = "ch1",
                                   true_betas = true_betas)
    
    # Run imputations, run cause-specific cox, and pool
    imp_ch12 <- mice(dat, m = m[length(m)],
                     method = method, 
                     predictorMatrix = mpred_ch12,
                     print = FALSE)
    
    results_ch12 <- mice_pool_diffm(imps = imp_ch12,
                                    m = m, j = j, label = "ch12",
                                    true_betas = true_betas)
    
    
    imp_ch12_int <- mice(dat, m = m[length(m)],
                         method = method,
                         predictorMatrix = mpred_ch12_int,
                         print = FALSE)
    
    results_ch12_int <- mice_pool_diffm(imps = imp_ch12_int, 
                                        m = m, j = j, label = "ch12_int",
                                        true_betas = true_betas)
    
    results_MI <- bind_rows(results_ch1, results_ch12, results_ch12_int) 
    
    return(results_MI) 
  })
  
  # Bring the rest together
  agg <- format_MI_repet(repetitions)
  
  
  # Bartlett smcfs model - removes warnings in case rej sampling fails
  mod_smcfcs <- suppressWarnings(
    quiet(smcfcs(originaldata = dat, 
                 smtype = "compet", 
                 smformula = c("Surv(t, eps == 1) ~ X1 + X2",
                               "Surv(t, eps == 2) ~ X1 + X2"), 
                 method = c("norm", rep("", 10)), 
                 m = m[length(m)]))
  )
    
  
  summ_smcfcs <- smcfcs_pool_diffm(mod_smcfcs, m, label = "smcfcs", 
                                   true_betas = true_betas) 
  
  ## Ref model: 
  mod_ref <- coxph(Surv(t, eps == 1) ~ X1_orig + X2, data = dat)
  
  # Read in summary, compute CI & label analysis type
  results_ref <- format_CCA_ref(mod_ref, true_betas = true_betas, label = "1ref")
  
  ## CCA:
  mod_CCA <- coxph(Surv(t, eps == 1) ~ X1 + X2, data = dat)
  
  results_CCA <- format_CCA_ref(mod_CCA, true_betas = true_betas, label = "CCA")
  
  # Results combined
  result <- rbind.data.frame(results_CCA, agg, summ_smcfcs) %>% 
    select(-pval, -`2.5 %`, -`97.5 %`)
                         
  # Label it 
  result$scen_name <- scenarios$scen_name[row]
  result$seed <- seed 
  
  return(result)
})

results <- bind_rows(res_seed)

#save(results, file = paste0("/exports/molepi/users/efbonneville/mi_comprisks/results_ISCB/results_", 
#                            seed, ".RData"))


# Uncomment and run below if editing file from gitlab directory - to go back to main
#setwd("../../code_MI_comprisks/")
#file.edit("simulations_main.R")
