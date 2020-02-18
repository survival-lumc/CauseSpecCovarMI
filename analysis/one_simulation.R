##******************************************##
## Test building for running one simulation ##
##******************************************##


# Data generation / parameter fixing --------------------------------------


# Fixed/varied parameters:

n <- 500 # sample size
meas <- "continous" # measurement level of X
rho <- 0.5 # correlation X and Z
beta1 <- 0.5 # Coeff of X in CauseSpec mod of event 1
mechan <- "MCAR" # missingness mechanism
prop_miss <- .05 # proportion of missingness
m <- c(2, 5) # Number of imputations of interest
iters_MI <- 5 # Iterations of multiple imputation procedure
horiz <- c(0.5, 1) # Prediction horizons
rate_cens <- 0.01 # Exponential censoring rate
eta1 <- NULL # Eta in MAR/MNAR/MAR GEN

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


# Set up imputations
mats <- get_predictor_mats(dat) # predictor matrix
methods <- get_imp_models(dat) # imputation models


# Ch1
imp_ch1 <- mice(dat, m = m[length(m)],
                method = methods, 
                predictorMatrix = mats$CH1,
                maxit = iters_MI, 
                print = FALSE)


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
imps_smcfcs <- quiet(
  record_warning(
    smcfcs(originaldata = dat, 
           smtype = "compet", 
           smformula = c("Surv(t, eps == 1) ~ X + Z",
                         "Surv(t, eps == 2) ~ X + Z"), 
           method = methods, 
           m = m[length(m)], 
           numit = iters_MI)
  )
)


# Extract rej sampling - add to summary data-frame after
str_extract(imps_smcfcs$warning, "[0-9]+")


# JointAI - make iterations to 1000
# Pick n.iter such that MCSE/SD smaller than 5% 
mod <- coxph_imp(Surv(t, ev1) ~ X + Z, 
                 n.iter = 25, 
                 data = dat, 
                 models = methods[methods != ""])

# 500 iters (500 adaptation) on datset of 500 takes: 4.65minutes
# [...] for size = 2000 : 

# Take mean and SD of posterior as estimate + SE
MC_error(mod)



# Summarise / pool estimates ----------------------------------------------


# Store all complete imputed dataset for mice/smcfcs
complist <- list("ch1" = mice::complete(imp_ch1, action = "all"),
                 "ch12" = mice::complete(imp_ch12, action = "all"),
                 "ch12_int" = mice::complete(imp_ch12_int, action = "all"),
                 "smcfcs" = imps_smcfcs$value$impDatasets)


# Run cox models on imputed datasets
mods_complist <- purrr::modify_depth(complist, .depth = 2, ~ setup_mstate(.x))  


# Pool estimates (mice, smcfcs, JointAI, CCA, ref)
estimates <- purrr::imap_dfr(mods_complist, 
                             ~ pool_diffm(.x, n_imp = m, analy = .y)) %>% 
  
  # Bind Bayes, CCA, ref
  bind_rows(
    summarise_ref_CCA(mod_ref, analy = "ref"),
    summarise_ref_CCA(mod_CCA, analy = "CCA"),
    summarise_bayes(mod)
  ) %>% 
  
  # Add true values
  mutate(true = case_when(
    str_detect(var, "X.1") ~ ev1_pars$b1,
    str_detect(var, "Z.1") ~ ev1_pars$gamm1,
    str_detect(var, "X.2") ~ ev2_pars$b2,
    str_detect(var, "Z.2") ~ ev2_pars$gamm2
  ))
  
View(estimates)


# Predicted probabilities -------------------------------------------------


# Make predictions for cox models fitted in each imputed dataset 
preds_list <- purrr::modify_depth(
  cox_mods, .depth = 2,
  ~ preds_mstate(cox_long = .x,
                 grid_obj = make_covar_grid(dat), 
                 times = horiz, 
                 ev1_pars = ev1_pars,
                 ev2_pars = ev2_pars)
)


# Pool predictions
pooled_preds <- purrr::imap_dfr(preds_list, 
                                ~ pool_diffm_preds(.x, n_imp = m, analy = .y))


View(pooled_preds)


# Scenario identifiers + seed + export to .Rdata --------------------------


# Scenario/seed should at least be in the filenames; also record smcfcs 
# Rej sampling failures too?





# Export to .Rdata


# -- End