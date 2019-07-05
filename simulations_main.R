#############################
## Main simulations script ##
#############################


# Clear environment
rm(list = ls())

# Load necessary packages
library(MASS)
library(tidyverse)
library(mice)
library(mitools)
library(survival)
library(smcfcs) # Bartlett package
#library(cmprsk)
#library(xtable)


# Read in support functions
source("data_generation/dat_generation_W.R")
source("support_functions/mice_pool_diffm.R")
source("support_functions/smcfcs_pool_diffm.R")
source("support_functions/quiet_printcat.R")


# Set a global seed
set.seed(1984)


# Fix parameters (equiv of scenarios later) :

N <- 2 # number of datasets to create
true_betas <- c(0.5, -0.5) # beta (X1) and gamma (X2) in data generating model
mech <- "MAR" # missingness mechanism
m <- c(1, 5, 10) # number of imputations to compare
method <- "norm" # imputation method, here bayesian linear regression
prob <- 0.3 # percentage missingness in X!
R <- 10 # Number times to replicate mice imputations


# Generate independent datasets & assign index
dats <- lapply(1:N, function(i) {
  
  dat <- dat_gener_W(N = 200, 
                     X_type = "contin",
                     mus = c(0, 0), 
                     covmat = matrix(c(1, 0.25, 
                                       0.25, 1), nrow = 2), # make into correlation mat
                     mech = mech, 
                     pars_MAR = c(1, 1),
                     p = prob, 
                     cause2 = "weib", 
                     vals_t1 = c(0.3, 1, true_betas), 
                     vals_t2 = c(1.7, 0.5, -0.5, 0.5)) 
  
  return(list("dat" = dat, "rep" = i))
})


# Assign MICE matrices once outside of lapply loop
dat_mats <- dats[[1]]$dat
mpred_ch1 <- mpred_ch12 <- mpred_ch12_int <-  matrix(0, ncol(dat_mats), ncol(dat_mats),
                                                     dimnames = list(names(dat_mats), names(dat_mats)))

## Ch1 model: MI with Z, eps (as factor) and H1(t) as covars
mpred_ch1["X1", c("X2", "ev1", "H1")] <- 1

## Ch12 model: MI with Z, eps, H1(t) & H2(t) as covars
mpred_ch12["X1", c("X2", "eps", "H1", "H2")] <- 1

# Ch12 model with addition of interactions H x Z
mpred_ch12_int["X1", c("X2", "eps", "H1", "H2", "H1_X2", "H2_X2")] <- 1

rm(dat_mats)


# Run simulations here
final_test <- lapply(dats, function(obj) {
  
  # Extract data set
  dat <- obj$dat
  
  # set.seed(obj$rep)
  
  # Start replicates for MI methods
  replicates <- lapply(1:R, function(r) {
    
    ## Multiple imputation section
    
    # Run imputations, run cause-specific cox, and pool - 
    # tidyversify this using https://stefvanbuuren.name/fimd/workflow.html
    imp_ch1 <- mice(dat, m = m[length(m)],
                    method = method, 
                    predictorMatrix = mpred_ch1,
                    print = FALSE)
    
    results_ch1 <- mice_pool_diffm(imps = imp_ch1,
                                   m = m, r = r, label = "ch1")
    
    # Run imputations, run cause-specific cox, and pool
    imp_ch12 <- mice(dat, m = m[length(m)],
                     method = method, 
                     predictorMatrix = mpred_ch12,
                     print = FALSE)
    
    results_ch12 <- mice_pool_diffm(imps = imp_ch12,
                                    m = m, r = r, label = "ch12")
    
    
    imp_ch12_int <- mice(dat, m = m[length(m)],
                         method = method,
                         predictorMatrix = mpred_ch12_int,
                         print = FALSE)
    
    results_ch12_int <- mice_pool_diffm(imps = imp_ch12_int, 
                                        m = m, r = r, label = "ch12_int")
    
    results_MI <- bind_rows(results_ch1, results_ch12, results_ch12_int) 
    
    return(results_MI) 
  })
  
  
  # keep SE from first replicate
  se_rep1 <- bind_rows(replicates) %>% 
    filter(rep == 1) %>% 
    select(se, var, analy, m)
  
  agg <- bind_rows(replicates) %>%
    
    # Power and coverage
    #mutate(pow = pval < 0.05,
    #       cov = `2.5 %` < coef && coef < `97.5 %`)
    group_by(var, analy, m) %>%
    mutate(sd_reps = sd(coef)) %>%
    summarise_all(~ round(mean(.), 3)) %>%
    ungroup() %>%
    select(-rep, -se) %>%
    right_join(se_rep1, by = c("var", "analy", "m")) %>%
    as.data.frame()
  
  
  # Bartlett smcfs model
  mod_smcfcs <- quiet(smcfcs(originaldata = dat, 
                             smtype = "compet", 
                             smformula = c("Surv(t, eps == 1) ~ X1 + X2",
                                           "Surv(t, eps == 2) ~ X1 + X2"), 
                             method = c("norm", rep("", 10)), 
                             m = m[length(m)]))
  
  summ_smcfcs <- smcfcs_pool_diffm(mod_smcfcs, m, label = "smcfcs") %>%
    mutate(sd_reps = 0)
                       
  
  ## Ref model: 
  mod_ref <- coxph(Surv(t, eps == 1) ~ X1_orig + X2, data = dat)
  
  # Read in summary, compute CI & label analysis type
  summ_ref <- summary(mod_ref)$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")]
  results_ref <- cbind.data.frame(summ_ref, confint(mod_ref)) %>%
    mutate(analy = "1ref")
  
  
  ## CCA:
  mod_CCA <- coxph(Surv(t, eps == 1) ~ X1 + X2, data = dat)
  
  # Read in summary
  summ_CCA <- summary(mod_CCA)$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")]
  
  # Combine / format results together with reference
  results_CCA <- cbind.data.frame(summ_CCA, confint(mod_CCA)) %>%
    mutate(analy = "CCA") %>% 
    bind_rows(results_ref) %>%
    mutate(var = rep(c("X1", "X2"), 2),
           true = rep(true_betas, 2),
           m = 0, pval = `Pr(>|z|)`,
           #sd_imps = 0,
           sd_reps = 0) %>%
    select(var, analy, m, coef, se = `se(coef)`, 
           pval, `2.5 %`, `97.5 %`, true, sd_reps) 
  
  
  # Make a counter
  #seq_perc <- seq(0, N, by = floor(N / 10))
  #if (obj$rep %in% seq_perc) {cat(100 * (obj$rep / N), "%")}
  
  return(rbind.data.frame(results_CCA, agg, summ_smcfcs))
})



results <- bind_rows(final_test) %>%
  group_by(var, analy, m) %>%
  mutate(se_emp = sd(coef)) %>%
  summarise_all(~ round(mean(.), 3)) %>%
  select(var, analy, m, coef, se, se_emp) %>%
  as.data.frame()


# For shark
#save(dat, file = "/exports/molepi/users/efbonneville/mi_comprisks/results.RData")


# After running 
# print(xtable(manip), digits = 3), include.rownames = FALSE)
