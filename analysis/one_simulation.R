##******************************************##
## Test building for running one simulation ##
##******************************************##


# Run system time on this whole thing

# Data generation / parameter fixing --------------------------------------


# Fixed/varied parameters:

n <- 500 # sample size
meas <- "binary" # measurement level of X
rho <- 0.5 # correlation X and Z
beta1 <- 1 # Coeff of X in CauseSpec mod of event 1
mechan <- "MNAR" # missingness mechanism
prop_miss <- .5 # proportion of missingness
m <- c(10) #25, 100) # Number of imputations of interest
iters_MI <- 50 # Iterations of multiple imputation procedure
horiz <- c(0.5, 1, 1.5) # Prediction horizons
rate_cens <- 0.01 # Exponential censoring rate
eta1 <- 2 # Eta in MAR/MNAR/MAR GEN


# Load scenarios RDS
scen_test <- scens %>% dplyr::slice(18)


one_simulation <- function(scenario, # scenario
                           rep_num) { # repetition number
  
  # Reproducibiliy
  seed <- scenario$seed + rep_num
  set.seed(seed)
  

  # Data generation section ----
  

  # Extract parameters from AFTs ran on MDS-long term data
  baseline <- readRDS(
    "analysis/data/derived_data/MDS_shape_rates.rds"
  )
  
  # Parameter Weibull event 1
  ev1_pars <- list(
    "a1" = baseline[baseline$state == "REL", "shape"], 
    "h1_0" = baseline[baseline$state == "REL", "rate"],
    "b1" = scenario$beta1, 
    "gamm1" = 1
  )
  
  # Parameters Weibull event 2
  ev2_pars <- list(
    "a2" = baseline[baseline$state == "NRM", "shape"], 
    "h2_0" = baseline[baseline$state == "NRM", "rate"], 
    "b2" = .5, 
    "gamm2" = .5
  )

  # Generate a dataset based on scenario
  dat <- generate_dat(n = scenario$n,
                      X_type = scenario$X_level, 
                      r = scenario$rho, 
                      ev1_pars = ev1_pars,
                      ev2_pars = ev2_pars, 
                      rate_cens = baseline[baseline$state == "EFS", "rate"], 
                      mech = scenario$miss_mech, 
                      p = scenario$prop_miss,
                      eta1 = scenario$eta1)
  
  
  # Run reference and CCA models ---- 
  
  
  # Reference (full dataset)
  mod_ref <- setup_mstate(dat %>% mutate(X = X_orig))
  
  # Complete case analyses
  mod_CCA <- setup_mstate(dat)
  
  
  # Imputation part ----
  
  
  # Fixed parameters
  m <- c(2, 3, 4) # Number of imputations of interest, should be c(2, 5, 25, 100)
  iters_MI <- 5 # Iterations of multiple imputation procedure, = 25
  horiz <- c(0.5, 2, 5) # Prediction horizons, 6mo, 2Y, 5Y
  
  # Set methods and predictor matrices
  mats <- get_predictor_mats(dat) 
  methods <- get_imp_models(dat) 
  
  # Ch1 model
  imp_ch1 <- mice(dat, m = m[length(m)],
                  method = methods, 
                  predictorMatrix = mats$CH1,
                  maxit = iters_MI, 
                  print = FALSE)
  
  # Ch12 model
  imp_ch12 <- mice(dat, m = m[length(m)],
                   method = methods, 
                   predictorMatrix = mats$CH12,
                   maxit = iters_MI, 
                   print = FALSE)
  
  # Ch12_int model
  imp_ch12_int <- mice(dat, m = m[length(m)],
                       method = methods, 
                       predictorMatrix = mats$CH12_int,
                       maxit = iters_MI, 
                       print = FALSE)
  
  
  # smcfcs 
  imp_smcfcs <- quiet(
    record_warning(
      smcfcs(originaldata = dat, 
             smtype = "compet", 
             smformula = c("Surv(t, eps == 1) ~ X + Z",
                           "Surv(t, eps == 2) ~ X + Z"), 
             method = methods, 
             m = m[length(m)], 
             numit = iters_MI, #iters_MI,
             rjlimit = 5000) # 5 times higher than default, avoid rej sampling errors
    )
  )
  
  # Extract rej sampling - add to summary data-frame after, columns called warns?
  str_extract(imp_smcfcs$warning, "[0-9]+")
  
  # Store all complete imputed datasets
  complist <- list("ch1" = mice::complete(imp_ch1, action = "all"),
                   "ch12" = mice::complete(imp_ch12, action = "all"),
                   "ch12_int" = mice::complete(imp_ch12_int, action = "all"),
                   "smcfcs" = imp_smcfcs$value$impDatasets)
  
  # Run cox models on imputed datasets
  mods_complist <- purrr::modify_depth(complist, .depth = 2, ~ setup_mstate(.x)) 
  
  
  # Summarise/pool regression coefficients ----
  
  
  estimates <- purrr::imap_dfr(mods_complist, 
                               ~ pool_diffm(.x, n_imp = m, analy = .y)) %>% 
    
    # Bind Bayes, CCA, ref
    bind_rows(
      summarise_ref_CCA(mod_ref, analy = "ref"),
      summarise_ref_CCA(mod_CCA, analy = "CCA")#,
      #summarise_bayes(mod)
    ) %>% 
    
    # Add true values
    mutate(true = case_when(
      str_detect(var, "X.1") ~ ev1_pars$b1,
      str_detect(var, "Z.1") ~ ev1_pars$gamm1,
      str_detect(var, "X.2") ~ ev2_pars$b2,
      str_detect(var, "Z.2") ~ ev2_pars$gamm2
    ))
  
  
  # Save as RDS...analysis/simulation results/estimates
  #saveRDS(estimates, file = "scenariolabel/number_repnum_estims")
  
  
  # Prediction part ----
  
  
  # Make predictions for cox models fitted in each imputed dataset 
  preds_list <- purrr::modify_depth(
    mods_complist, .depth = 2,
    ~ preds_mstate(cox_long = .x,
                   grid_obj = make_covar_grid(dat), 
                   times = horiz, 
                   ev1_pars = ev1_pars,
                   ev2_pars = ev2_pars)
  )
  
  
  # Pool predictions
  pooled_preds <- purrr::imap_dfr(preds_list, 
                                  ~ pool_diffm_preds(.x, n_imp = m, analy = .y))
  
  # Save as RDS in analysis/simulation results/predictions
  #saveRDS(pooled_preds, file = "scenariolabel/number_repnum_preds")
  
  # Also store seeds
  
  # Remove this return thing after
  return(list("estims" = estimates,
              "preds" = pooled_preds))
}


test_run <- one_simulation(scen_test, rep_num = 1)

View(test_run$estims %>%  dplyr::mutate_if(is.numeric, ~ round(., 3)))



# Fixed parameters


# Parameter Weibull event 1
ev1_pars = list("a1" = 2, "h1_0" = 1,
                "b1" = beta1, "gamm1" = 1)

# Parameters Weibull event 2
ev2_pars = list("a2" = 2.5, "h2_0" = .5, 
                "b2" = .5, "gamm2" = .5)

# Set some seed - depending on scenario / repetition
set.seed(1)

# Generate a dataset
dat <- generate_dat(n = n,
                    X_type = meas, 
                    r = rho, 
                    ev1_pars = ev1_pars,
                    ev2_pars = ev2_pars, 
                    rate_cens = rate_cens, # fixed from data
                    mech = mechan, 
                    p = prop_miss,
                    eta1 = eta1)

# Run all methods ---------------------------------------------------------


# Reference (full dataset)
mod_ref <- setup_mstate(dat %>% mutate(X = X_orig))


# Complete case analyses
mod_CCA <- setup_mstate(dat)

# Small plot for fun
#cumincs_plot_truepred(mod_CCA, 
#                      combo = data.frame("val_X" = 0,
#                                        "val_Z" = 0), 
#                     ev1_pars = ev1_pars,
#                      ev2_pars = ev2_pars, 
#                      dat = dat)


# Set up imputations
mats <- get_predictor_mats(dat) # predictor matrix
methods <- get_imp_models(dat) # imputation models


# Ch1 - about 20 seconds for nimp = 100, and iters = 25
imp_ch1 <- mice(dat, m = m[length(m)],
                method = methods, 
                predictorMatrix = mats$CH1,
                maxit = iters_MI, 
                print = FALSE)


purrr::map(mice::complete(imp_ch1, action = "all"), 
           ~ setup_mstate(.x)) %$%
  pool_diffm(., m, analy = "ch1")


# Ch12
imp_ch12 <- mice(dat, m = m[length(m)],
                 method = methods, 
                 predictorMatrix = mats$CH12,
                 maxit = iters_MI, 
                 print = FALSE)

# Ch12_int
imp_ch12_int <- mice(dat, m = m[length(m)],
                     method = methods, 
                     predictorMatrix = mats$CH12_int,
                     maxit = iters_MI, 
                     print = FALSE)


# smcfcs - wrap with quiet() after, capture rejection sampling errors
# Check if methods ok for log reg
# Takes about 15.5 minutes seconds for nimp = 100, and iters = 25 
imp_smcfcs <- quiet(
  record_warning(
    smcfcs(originaldata = dat, 
           smtype = "compet", 
           smformula = c("Surv(t, eps == 1) ~ X + Z",
                         "Surv(t, eps == 2) ~ X + Z"), 
           method = methods, 
           m = m[length(m)], 
           numit = iters_MI, #iters_MI,
           rjlimit = 5000) # 5 times higher than default, avoid rej sampling errors
  )
)


# Extract rej sampling - add to summary data-frame after
str_extract(imp_smcfcs$warning, "[0-9]+")

purrr::map_dfr(1:m[length(m)], function(i) {
  as.data.frame(t(imp_smcfcs$value$smCoefIter[i, ,])) %>% 
    mutate(iter = 1:iters_MI)
}, .id = "imp_dat") %>% 
  rename("X1" = V1, "X2" = V2,
         "Z1" = V3, "Z2" = V4) %>% 
  gather("var", "value", X1:Z2) %>% 
  ggplot(aes(iter, value, col = imp_dat)) +
  geom_line() +
  theme(legend.position = "none") +
  facet_wrap(~ var)




# JointAI - make iterations to 1000
# Pick n.iter such that MCSE/SD smaller than 5% 
# This is not correct -
# See Bartlett comp risks paper - bayesian equiv of FCS survival

#mod <- coxph_imp(Surv(t, ev1) ~ X + Z, 
#                 n.iter = 1000, 
#                data = dat, 
#                 models = methods[methods != ""])

# 500 iters (500 adaptation) on datset of 500 takes: 4.65minutes
# [...] for size = 2000, iters = 1000 takes: 33.44 minutes :o
# with MCSE/SD or around 0.04-0.05

# Take mean and SD of posterior as estimate + SE
#MC_error(mod)



# Summarise / pool estimates ----------------------------------------------


# Store all complete imputed dataset for mice/smcfcs
complist <- list("ch1" = mice::complete(imp_ch1, action = "all"),
                 "ch12" = mice::complete(imp_ch12, action = "all"),
                 "ch12_int" = mice::complete(imp_ch12_int, action = "all"),
                 "smcfcs" = imp_smcfcs$value$impDatasets)


# Run cox models on imputed datasets
mods_complist <- purrr::modify_depth(complist, .depth = 2, ~ setup_mstate(.x))  


# Pool estimates (mice, smcfcs, JointAI, CCA, ref)
estimates <- purrr::imap_dfr(mods_complist, 
                             ~ pool_diffm(.x, n_imp = m, analy = .y)) %>% 
  
  # Bind Bayes, CCA, ref
  bind_rows(
    summarise_ref_CCA(mod_ref, analy = "ref"),
    summarise_ref_CCA(mod_CCA, analy = "CCA")#,
    #summarise_bayes(mod)
  ) %>% 
  
  # Add true values
  mutate(true = case_when(
    str_detect(var, "X.1") ~ ev1_pars$b1,
    str_detect(var, "Z.1") ~ ev1_pars$gamm1,
    str_detect(var, "X.2") ~ ev2_pars$b2,
    str_detect(var, "Z.2") ~ ev2_pars$gamm2
  ))
  
View(estimates %>% 
       mutate_if(is.numeric, ~ round(., 3)))


# Predicted probabilities -------------------------------------------------


# Make predictions for cox models fitted in each imputed dataset 
preds_list <- purrr::modify_depth(
  mods_complist, .depth = 2,
  ~ preds_mstate(cox_long = .x,
                 grid_obj = make_covar_grid(dat), 
                 times = horiz, 
                 ev1_pars = ev1_pars,
                 ev2_pars = ev2_pars)
)


# Pool predictions
pooled_preds <- purrr::imap_dfr(preds_list, 
                                ~ pool_diffm_preds(.x, n_imp = m, analy = .y))


View(pooled_preds %>% 
       mutate_if(is.numeric, ~ round(., 3)))


# Scenario identifiers + seed + export to .Rdata --------------------------


# Scenario/seed should at least be in the filenames; also record smcfcs 
# Rej sampling failures too?





# Export to .Rdata


# -- End