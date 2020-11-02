##*************************************##
## One replication of one sim scenario ##
##      for ordered categorical X      ##
##*************************************##

# This will go to run simulation
devtools::load_all()
scenarios <- readRDS("./inst/testdata/scenarios_ordcat.rds") 

scenario <- scenarios[scenarios$scen_num == 134, ]
rep_num <- 1

one_simulation_ordcat <- function(scenario, rep_num) {
  
  # Reproducibiliy
  seed <- scenario$seed + rep_num
  set.seed(seed)
  
  # Data generation section ----
  
  
  # Extract parameters from AFTs ran on MDS-long term data
  baseline <- readRDS(
    "./inst/testdata/MDS_shape_rates.rds"
  )
  
  # Check shape of 1st hazard
  if (scenario$haz_shape == "similar") {
    shape_ev1 <- baseline[baseline$state == "REL", "shape"]
    base_rate_ev1 <- baseline[baseline$state == "REL", "rate"]
    
  } else {
    
    # shape > 1 so increasing hazard -> lower base rate to 0.04
    # which avoids 10y EFS of zero!
    shape_ev1 <- 1.5
    base_rate_ev1 <- 0.04
  }
  
  # Parameter Weibull event 1
  ev1_pars <- list(
    "a1" = shape_ev1, 
    "h1_0" = base_rate_ev1,
    "b1" = c(0.75, scenario$beta1), 
    "gamm1" = 1
  )
  
  # Parameters Weibull event 2
  ev2_pars <- list(
    "a2" = baseline[baseline$state == "NRM", "shape"], 
    "h2_0" = baseline[baseline$state == "NRM", "rate"], 
    "b2" = c(0.5, 0.75), 
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
  
  # Delete after
  dat %>% 
    ggplot(aes(X_orig, Z, col = factor(miss_ind))) +
    geom_jitter(width = 0.2, alpha = 0.75)
  
  # Run reference and CCA models ---- 
  
  
  # Reference (full dataset)
  mod_ref <- SimsCauseSpecCovarMiss::setup_mstate(
    dat %>% dplyr::mutate(X = X_orig)
  )
  
  # Complete case analyses
  mod_CCA <- SimsCauseSpecCovarMiss::setup_mstate(dat)
  
  
  # Imputation part ----
  
  
  # Fixed parameters
  m <- c(2, 3, 4) #c(5, 10, 25, 50) # Number of imputations of interest
  iters_MI <- 3 #20 # Iterations of multiple imputation procedure
  
  # Set methods and predictor matrices
  mats <- SimsCauseSpecCovarMiss::get_predictor_mats(dat) 
  
  # For binary mice = smcfcs
  methods_mice <- SimsCauseSpecCovarMiss::set_mi_methods(
    dat = dat, 
    var_names_miss = naniar::miss_var_which(dat), 
    imp_type = "mice", 
  ) 
  
  methods_smcfcs <- SimsCauseSpecCovarMiss::set_mi_methods(
    dat = dat, 
    var_names_miss = naniar::miss_var_which(dat), 
    imp_type = "smcfcs", 
  ) 
  
  cat("Running imp_ch1... \n\n")
  
  # Ch1 model
  imp_ch1 <- mice::mice(dat, m = m[length(m)],
                        method = methods_mice, 
                        predictorMatrix = mats$CH1,
                        maxit = iters_MI, 
                        print = FALSE)
  
  cat("\n Running imp_ch12... \n\n")
  
  
  # Ch12 model
  imp_ch12 <- mice::mice(dat, m = m[length(m)],
                         method = methods_mice, 
                         predictorMatrix = mats$CH12,
                         maxit = iters_MI, 
                         print = FALSE, 
                         threshold = 1) # Avoid removing H2 
  
  cat("\n Running imp_ch12_int... \n\n")
  
  # Ch12_int model
  imp_ch12_int <- mice::mice(dat, m = m[length(m)],
                             method = methods_mice, 
                             predictorMatrix = mats$CH12_int,
                             maxit = iters_MI, 
                             print = FALSE, 
                             threshold = 1) # Avoid removing H2, and H2_Z 
  
  cat("\n Running smcfcs... \n\n")
  
  # Smcfcs - quiet stops printing
  imp_smcfcs <- #SimsCauseSpecCovarMiss::quiet(
    
    # Record number of rej sampling failures, if any
    SimsCauseSpecCovarMiss::record_warning(
      
      # Smcfcs starts here
      smcfcs::smcfcs(
        originaldata = dat, 
        smtype = "compet", 
        smformula = c("Surv(t, eps == 1) ~ X + Z",
                      "Surv(t, eps == 2) ~ X + Z"), 
        method = methods_smcfcs, 
        m = m[length(m)], 
        numit = iters_MI, 
        rjlimit = 5000) # 5 times higher than default, avoid rej sampling errors
    )
  #)
  
  # Store all complete imputed datasets
  complist <- list("ch1" = mice::complete(imp_ch1, action = "all"),
                   "ch12" = mice::complete(imp_ch12, action = "all"),
                   "ch12_int" = mice::complete(imp_ch12_int, action = "all"),
                   "smcfcs" = imp_smcfcs$value$impDatasets)
  
  cat("\n Fitting mstate model in all 200 imputed datasets... \n\n")
  
  # Run cox models on imputed datasets
  mods_complist <- purrr::modify_depth(complist, .depth = 2, ~ setup_mstate(.x)) 
  
  
  # Summarise/pool regression coefficients ----
  
  cat("\n Pooling estimates... \n\n")
  
  # Make this into a function?
  estimates <- purrr::imap_dfr(
    mods_complist, 
    ~ SimsCauseSpecCovarMiss::pool_diffm(.x, n_imp = m, analy = .y)
  ) %>% 
    
    # Bind Bayes, CCA, ref
    dplyr::bind_rows(
      SimsCauseSpecCovarMiss::summarise_ref_CCA(mod_ref, analy = "ref"),
      SimsCauseSpecCovarMiss::summarise_ref_CCA(mod_CCA, analy = "CCA")
    ) %>% 
    
    # Add true values
    dplyr::mutate(true = dplyr::case_when(
      stringr::str_detect(var, "Xinterm.1") ~ ev1_pars$b1[1],
      stringr::str_detect(var, "Xhigh.1") ~ ev1_pars$b1[2],
      stringr::str_detect(var, "Z.1") ~ ev1_pars$gamm1,
      stringr::str_detect(var, "Xinterm.2") ~ ev2_pars$b2[1],
      stringr::str_detect(var, "Xhigh.2") ~ ev2_pars$b2[2],
      stringr::str_detect(var, "Z.2") ~ ev2_pars$gamm2
    )) %>% 
    
    # Add mice warnings and rej sampling errors for smcfcs
    dplyr::mutate(
      warns = dplyr::case_when(
        analy == "smcfcs" ~ as.numeric(stringr::str_extract(imp_smcfcs$warning,
                                                            "[0-9]+")),
        analy == "ch12" ~ as.numeric(!(is.null(imp_ch12$loggedEvents))),
        analy == "ch12_int" ~ as.numeric(!(is.null(imp_ch12_int$loggedEvents))),
        analy %in% c("CCA", "ref", "ch1") ~ 0
      ) 
    ) %>% 
    
    # Add scenario summary
    cbind.data.frame(scen_summary = add_scen_details(
      scenario, 
      seed = seed, 
      rep_num = rep_num
    ), row.names = NULL, 
    stringsAsFactors = FALSE)
  
  # RDS saving and storing seeds/scen_num ----
  
  # Save as RDS in analysis/simulation results/predictions
  # saveRDS(
  #   estimates, 
  #   file = paste0("./analysis/simulations/sim-reps_individual/estims_scen",
  #                 scenario$scen_num, "_seed", seed, ".rds")
  # )
  
  # No prediction part - do to too many possible ref patients
  
}