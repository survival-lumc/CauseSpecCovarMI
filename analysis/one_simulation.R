##******************************************##
## Test building for running one simulation ##
##******************************************##


# Run system time on this whole thing

# Load scenarios RDS
scenario <- scens %>% dplyr::slice(4)
scenario
rep_num <- 1

one_simulation <- function(scenario, # scenario
                           rep_num) { # repetition number
  
  # Reproducibiliy
  seed <- scenario$seed + rep_num
  set.seed(seed)
  

  # Data generation section ----
  

  # Extract parameters from AFTs ran on MDS-long term data
  baseline <- readRDS(
    system.file("testdata", 
                "MDS_shape_rates.rds", 
                package = "SimsCauseSpecCovarMiss")
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
  dat <- SimsCauseSpecCovarMiss::generate_dat(
    n = scenario$n,
    X_type = scenario$X_level, 
    r = scenario$rho, 
    ev1_pars = ev1_pars,
    ev2_pars = ev2_pars, 
    rate_cens = baseline[baseline$state == "EFS", "rate"], 
    mech = scenario$miss_mech, 
    p = scenario$prop_miss,
    eta1 = scenario$eta1
  )
  
  
  # Run reference and CCA models ---- 
  
  
  # Reference (full dataset)
  mod_ref <- SimsCauseSpecCovarMiss::setup_mstate(
    dat %>% dplyr::mutate(X = X_orig)
  )
  
  # Complete case analyses
  mod_CCA <- SimsCauseSpecCovarMiss::setup_mstate(dat)
  
  
  # Imputation part ----
  
  
  # Fixed parameters
  m <- c(2, 3, 4) # Number of imputations of interest, should be c(2, 5, 25, 100)
  iters_MI <- 5 # Iterations of multiple imputation procedure, = 25
  
  # Set methods and predictor matrices
  mats <- SimsCauseSpecCovarMiss::get_predictor_mats(dat) 
  methods <- SimsCauseSpecCovarMiss::get_imp_models(dat) 
  
  # Ch1 model
  imp_ch1 <- mice::mice(dat, m = m[length(m)],
                        method = methods, 
                        predictorMatrix = mats$CH1,
                        maxit = iters_MI, 
                        print = FALSE)
  
  # Ch12 model
  imp_ch12 <- mice::mice(dat, m = m[length(m)],
                         method = methods, 
                         predictorMatrix = mats$CH12,
                         maxit = iters_MI, 
                         print = FALSE)
                   
  
  # Ch12_int model
  imp_ch12_int <- mice::mice(dat, m = m[length(m)],
                             method = methods, 
                             predictorMatrix = mats$CH12_int,
                             maxit = iters_MI, 
                             print = FALSE)
  
  # Smcfcs - quiet stops printing
  imp_smcfcs <- SimsCauseSpecCovarMiss::quiet(
    
    # Record number of rej sampling failures, if any
    SimsCauseSpecCovarMiss::record_warning(
      
      # Smcfcs starts here
      SimsCauseSpecCovarMiss::smcfcs_timefix(
        originaldata = dat, 
        smtype = "compet", 
        smformula = c("Surv(t, eps == 1) ~ X + Z",
                      "Surv(t, eps == 2) ~ X + Z"), 
        method = methods, 
        m = m[length(m)], 
        numit = iters_MI, 
        rjlimit = 5000) # 5 times higher than default, avoid rej sampling errors
    )
  )
  
  # Store all complete imputed datasets
  complist <- list("ch1" = mice::complete(imp_ch1, action = "all"),
                   "ch12" = mice::complete(imp_ch12, action = "all"),
                   "ch12_int" = mice::complete(imp_ch12_int, action = "all"),
                   "smcfcs" = imp_smcfcs$value$impDatasets)
  
  # Run cox models on imputed datasets
  mods_complist <- purrr::modify_depth(complist, .depth = 2, ~ setup_mstate(.x)) 
  
  
  # Summarise/pool regression coefficients ----
  
  
  estimates <- purrr::imap_dfr(
    mods_complist, 
    ~ SimsCauseSpecCovarMiss::pool_diffm(.x, n_imp = m, analy = .y)
  ) %>% 
    
    # Bind Bayes, CCA, ref
    dplyr::bind_rows(summarise_ref_CCA(mod_ref, analy = "ref"),
                     summarise_ref_CCA(mod_CCA, analy = "CCA")) %>% 
    
    # Add true values
    dplyr::mutate(true = dplyr::case_when(
      stringr::str_detect(var, "X.1") ~ ev1_pars$b1,
      stringr::str_detect(var, "Z.1") ~ ev1_pars$gamm1,
      stringr::str_detect(var, "X.2") ~ ev2_pars$b2,
      stringr::str_detect(var, "Z.2") ~ ev2_pars$gamm2
    )) %>% 
    
    # Add rej sampling errors for smcfcs
    dplyr::mutate(
      warnings = ifelse(
        analy == "smcfcs",
        as.numeric(str_extract(imp_smcfcs$warning, "[0-9]+")),
        0
      )
    )
  
  
  # Prediction part ----
  
  
  horiz <- c(0.5, 5, 10) # Prediction horizons, 6mo, 5Y, 10Y
  covar_grid <- SimsCauseSpecCovarMiss::make_covar_grid(dat)
  
  # Make predictions for cox models fitted in each imputed dataset 
  preds_list <- purrr::modify_depth(
    mods_complist, .depth = 2,
    ~ SimsCauseSpecCovarMiss::get_preds_grid(
        cox_long = .x,
        grid_obj = covar_grid, 
        times = horiz, 
        ev1_pars = ev1_pars,
        ev2_pars = ev2_pars
    )
  )
  
  
  # Pool predictions
  pooled_preds <- purrr::imap_dfr(
    preds_list, 
    ~ SimsCauseSpecCovarMiss::pool_diffm_preds(.x, n_imp = m, analy = .y)
  )
  

  # RDS saving and storing seeds/scen_num ----

  
  # Scenario/seed should at least be in the filenames
  
  # Save as RDS in analysis/simulation results/predictions
  #saveRDS(pooled_preds, file = "scenariolabel/number_repnum_preds")
  
  # For shark
  #saveRDS(dat, 
  # file = "/exports/molepi/users/efbonneville/.../filename/rds")
  
  
  # Remove this return thing after
  return(list("estims" = estimates,
              "preds" = pooled_preds))
}


# Do the timings!!
test_run <- one_simulation(scen_test, rep_num = 1)

View(test_run$estims %>%  dplyr::mutate_if(is.numeric, ~ round(., 3)))







# -- End

