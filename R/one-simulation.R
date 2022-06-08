##*************************************##
## One replication of one sim scenario ##
##*************************************##


#' Run one replication of one scenario
#' 
#' @param scenario A row from the scenarios dataframe (link)
#' @param rep_num Scalar, replication numbr
#' 
#' @export
one_simulation <- function(scenario, # scenario
                           rep_num) { # repetition number
  
  # Set seed based on scenario and repetition number
  seed <- scenario$seed + rep_num
  set.seed(seed)
  

  # Data generation section ----
  

  # Read parameters from AFTs ran on MDS-long term data
  baseline <- CauseSpecCovarMI::mds_shape_rates
  
  # Check shape of 1st hazard
  if (scenario$haz_shape == "similar") {
    shape_ev1 <- baseline[baseline$state == "REL", "shape"]
    base_rate_ev1 <- baseline[baseline$state == "REL", "rate"]
    
  } else {
    
    # shape > 1 so increasing hazard -> lower base rate to 0.04 
    # it avoids 10y EFS of zero
    shape_ev1 <- 1.5
    base_rate_ev1 <- 0.04
  }
  
  # Parameter Weibull event 1
  ev1_pars <- list(
    "a1" = shape_ev1, 
    "h1_0" = base_rate_ev1,
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
  dat <- generate_dat(
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
  
  
  # Reference model (full dataset)
  mod_ref <- setup_mstate(dat %>% dplyr::mutate(X = .data$X_orig))
  
  # Complete case analyses
  mod_CCA <- setup_mstate(dat)
  
  
  # Imputation part ----
  
  
  # Fixed parameters
  m <- c(5, 10, 25, 50) # Number of imputations of interest
  iters_MI <- 1 #20 # Iterations of multiple imputation procedure
  # Can be set to 1 since chained equations are not needed, but original simulations
  # were unfortunately run with 20 
  
  # Set methods and predictor matrices
  mats <- get_predictor_mats(dat) 
  
  # For binary mice = smcfcs
  methods <- set_mi_methods(
    dat = dat, 
    var_names_miss = naniar::miss_var_which(dat), 
    imp_type = "mice", 
    cont_method = "norm" 
  ) 
  
  cat("Running imp_ch1... \n\n")
  
  # Ch1 model
  imp_ch1 <- mice::mice(
    dat, 
    m = m[length(m)],
    method = methods, 
    predictorMatrix = mats$CH1,
    maxit = iters_MI, 
    print = FALSE
  )
  
  cat("\n Running imp_ch12... \n\n")
  
  
  # Ch12 model
  imp_ch12 <- mice::mice(
    dat, 
    m = m[length(m)],
    method = methods, 
    predictorMatrix = mats$CH12,
    maxit = iters_MI, 
    print = FALSE, 
    threshold = 1
  ) # Avoid removing H2 
                   
  cat("\n Running imp_ch12_int... \n\n")
  
  # Ch12_int model
  imp_ch12_int <- mice::mice(
    dat, 
    m = m[length(m)],
    method = methods, 
    predictorMatrix = mats$CH12_int,
    maxit = iters_MI, 
    print = FALSE, 
    threshold = 1
  ) # Avoid removing H2, and H2_Z 
  
  cat("\n Running smcfcs... \n\n")
  
  # Smcfcs - 
  imp_smcfcs <- record_warning(
    smcfcs::smcfcs(
      originaldata = dat, 
      smtype = "compet", 
      smformula = c("Surv(t, eps == 1) ~ X + Z",
                    "Surv(t, eps == 2) ~ X + Z"), 
      method = methods, 
      m = m[length(m)], 
      numit = iters_MI, 
      rjlimit = 5000
    ) # 5 times higher than default, avoid rej sampling errors
  )
  
  # Store all complete imputed datasets
  complist <- list(
    "ch1" = mice::complete(imp_ch1, action = "all"),
    "ch12" = mice::complete(imp_ch12, action = "all"),
    "ch12_int" = mice::complete(imp_ch12_int, action = "all"),
    "smcfcs" = imp_smcfcs$value$impDatasets
  )
  
  cat("\n Fitting mstate model in all 200 imputed datasets... \n\n")
  
  # Run cox models on imputed datasets
  mods_complist <- purrr::modify_depth(complist, .depth = 2, ~ setup_mstate(.x)) 
  
  
  # Summarise/pool regression coefficients ----
  
  cat("\n Pooling estimates... \n\n")
  
  # Bind all estimates
  estimates <- purrr::imap_dfr(
    mods_complist, 
    ~ pool_diffm(.x, n_imp = m, analy = .y)
  ) %>% 
    
    # Bind CCA, ref
    dplyr::bind_rows(
      summarise_ref_CCA(mod_ref, analy = "ref"),
      summarise_ref_CCA(mod_CCA, analy = "CCA")
    ) %>% 
    
    # Add true values
    dplyr::mutate(
      true = dplyr::case_when(
        stringr::str_detect(var, "X.1") ~ ev1_pars$b1,
        stringr::str_detect(var, "Z.1") ~ ev1_pars$gamm1,
        stringr::str_detect(var, "X.2") ~ ev2_pars$b2,
        stringr::str_detect(var, "Z.2") ~ ev2_pars$gamm2
      )
    ) %>% 
    
    # Add mice warnings and rej sampling errors for smcfcs
    dplyr::mutate(
      warns = dplyr::case_when(
        analy == "smcfcs" ~ as.numeric(stringr::str_extract(imp_smcfcs$warning, "[0-9]+")),
        analy == "ch12" ~ as.numeric(!is.null(imp_ch12$loggedEvents)),
        analy == "ch12_int" ~ as.numeric(!(is.null(imp_ch12_int$loggedEvents))),
        analy %in% c("CCA", "ref", "ch1") ~ 0
      ) 
    ) %>% 
    
    # Add scenario summary
    cbind.data.frame(
      scen_summary = add_scen_details(
        scenario = scenario, 
        seed = seed, 
        rep_num = rep_num
      ), 
      row.names = NULL, 
      stringsAsFactors = FALSE
    )
  
  
  # Prediction part ----
  
  
  horiz <- c(0.5, 5, 10) # Prediction horizons, 6mo, 5Y, 10Y, c(0.5, 5, 10)
  covar_grid <- make_covar_grid(dat)
  
  # Compute the true ones once here to solve bottleneck?
  true_cis <- purrr::map(seq_len(nrow(covar_grid)), .f = function(row) {
    get_true_cuminc(ev1_pars, ev2_pars, times = horiz, combo = covar_grid[row, ])
  })
  
  cat("\n Making predictions and pooling... \n\n")
  
  # Make predictions for cox models fitted in each imputed dataset 
  preds_list <- purrr::modify_depth(
    mods_complist, .depth = 2,
    ~ {
      get_preds_grid(
        cox_long = .x,
        grid_obj = covar_grid, 
        times = horiz, 
        ev1_pars = ev1_pars,
        ev2_pars = ev2_pars,
        true_cuminc = true_cis
      )
    } 
  )
  
  # Make predictions also for ref and CCA
  preds_CCA <- get_preds_grid(
    cox_long = mod_CCA,
    grid_obj = covar_grid, 
    times = horiz, 
    ev1_pars = ev1_pars,
    ev2_pars = ev2_pars,
    true_cuminc = true_cis
  )
  
  preds_ref <- get_preds_grid(
    cox_long = mod_ref,
    grid_obj = covar_grid, 
    times = horiz, 
    ev1_pars = ev1_pars,
    ev2_pars = ev2_pars,
    true_cuminc = true_cis
  )
  
  
  # Pool predictions
  pooled_preds <- purrr::imap_dfr(
    preds_list, 
    ~ pool_diffm_preds(.x, n_imp = m, analy = .y)
  ) %>% 
    
    # Bind preds for CCA and ref
    dplyr::bind_rows(
      preds_CCA_ref(preds_CCA, "CCA"),
      preds_CCA_ref(preds_ref, "ref")
    ) %>% 
    
    # Add scenario summary label
    cbind.data.frame(
      scen_summary = add_scen_details(
        scenario, 
        seed = seed, 
        rep_num = rep_num
      ), 
      row.names = NULL, 
      stringsAsFactors = FALSE
    )
  

  cat("\n Replication done!")
  
  return(list("ests" = estimates, "preds" = pooled_preds))
  
  # RDS saving and storing seeds/scen_num ----
  
  # Save as RDS in analysis/simulation results/predictions
  # saveRDS( 
  #   object = estimates, 
  #   file = paste0(
  #     "data/sim-reps_indiv/regr/regr_scen", 
  #     scenario$scen_num, "_seed", seed, ".rds"
  #   )
  # )
  # saveRDS(
  #   object = pooled_preds, 
  #   file = paste0(
  #     "data/sim-reps_indiv/preds/preds_scen", 
  #     scenario$scen_num, "_seed", seed, ".rds"
  #   )
  # )
}

