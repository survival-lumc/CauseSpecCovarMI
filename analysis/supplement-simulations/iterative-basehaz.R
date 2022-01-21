# References and load libraries -------------------------------------------

devtools::load_all()
library(mice)
library(survival)
source("analysis/supplement-simulations/scenarios-supplemental.R")

# https://github.com/lbeesleyBIOSTAT/SRMIMI_Example_Code/blob/main/Example_Code_Normal.R
# https://github.com/alexanderrobitzsch/miceadds
# https://www.gerkovink.com/miceVignettes/Passive_Post_processing/Passive_imputation_post_processing.html

# One sim func ------------------------------------------------------------

# Make function to update basehaz
update_basehaz <- function(time, delta, x, z) {
  mod <- survival::coxph(
    Surv(time, delta) ~ x + z, 
    control = survival::coxph.control(timefix = FALSE)
  )
  basehaz_df <- survival::basehaz(mod, centered = FALSE)
  haz <- basehaz_df[match(time, basehaz_df[["time"]]), ][["hazard"]]
  return(haz)
}

one_simulation_breslow <- function(scenario, n_cpu, rep_num) {
  
  # Set the seed
  scen_number <- as.numeric(rownames(scenario))
  set.seed(-rep_num * scen_number)
  
  # Set up parameters
  baseline <- CauseSpecCovarMI::mds_shape_rates
  
  # Different basehaz setting
  alph1 <- 1.5; lam1 <- 0.04
  ev1_pars <- list("a1" = alph1, "h1_0" = lam1, "b1" = scenario$beta1, "gamm1" = 1)
  
  ev2_pars <- list(
    "a2" = baseline[baseline$state == "NRM", "shape"], 
    "h2_0" = baseline[baseline$state == "NRM", "rate"], 
    "b2" = .5, 
    "gamm2" = .5
  )
  
  # Generate dataset
  dat <- generate_dat(
    n = scenario$n,
    X_type = scenario$X_level, 
    r = 0.5, 
    ev1_pars = ev1_pars,
    ev2_pars = ev2_pars, 
    rate_cens = baseline[baseline$state == "EFS", "rate"], 
    mech = scenario$miss_mech, 
    p = scenario$prop_miss,
    eta1 = scenario$eta1
  )
  
  # Add true baseline hazards, and interactions
  dat$H1_true <- with(dat,  lam1 * t^alph1)
  dat$H2_true <- with(
    dat,
    baseline[baseline$state == "NRM", "rate"] * 
      t^(baseline[baseline$state == "NRM", "shape"])
  ) 
  dat$H1_true_Z <- with(dat, H1_true * Z)
  dat$H2_true_Z <- with(dat, H2_true * Z)
  
  # -- Run full dataset and CCA
  
  mod_ref <- setup_mstate(dat %>% dplyr::mutate(X = .data$X_orig))
  mod_CCA <- setup_mstate(dat)
  
  # -- smcfcs
  
  m <- 50 # make 25 or 50
  meths_smcfcs <- mice::make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
  
  imp_smcfcs <- record_warning(
    smcfcs::smcfcs.parallel(
      seed = round(runif(1, 0, 1e5)),
      n_cores = n_cpu,
      cl_type = "FORK",
      originaldata = dat, 
      smtype = "compet", 
      smformula = c(
        "Surv(t, eps == 1) ~ X + Z",
        "Surv(t, eps == 2) ~ X + Z"
      ), 
      method = meths_smcfcs, 
      m = m, 
      numit = 20, 
      rjlimit = 5000
    ) # 5 times higher than default, avoid rej sampling errors
  )
  
  #plot(imp_smcfcs$value)
  
  # -- mice setup
  
  mat <- mice::make.predictorMatrix(dat) 
  mat[] <- 0L 
  mat_nels <- mat_true <- mat_breslow <- mat 
  
  mat_nels["X", c("Z", "eps", "H1", "H2")] <- 1
  mat_true["X", c("Z", "eps", "H1_true", "H2_true")] <- 1
  meths_mice <- mice::make.method(dat, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
  
  # -- mice with nels and true
  
  imp_nels <- mice::mice(
    dat, 
    m = m,
    method = meths_mice, 
    predictorMatrix = mat_nels,
    maxit = 1, 
    print = FALSE
  )
  
  imp_true <- mice::mice(
    dat, 
    m = m,
    method = meths_mice, 
    predictorMatrix = mat_true,
    maxit = 1, 
    print = FALSE
  )
  
  # -- mice nels and true with interactions
  
  mat_nels["X", c("H1_Z", "H2_Z")] <- 1
  
  imp_nels_int <- mice::mice(
    dat, 
    m = m,
    method = meths_mice, 
    predictorMatrix = mat_nels,
    maxit = 1, 
    print = FALSE
  )
  
  mat_true["X", c("H1_true_Z", "H2_true_Z")] <- 1
  
  imp_true_int <- mice::mice(
    dat, 
    m = m,
    method = meths_mice, 
    predictorMatrix = mat_true,
    maxit = 1, 
    print = FALSE
  )
  
  # -- mice iterative breslow
  
  cat("\n Start of iterative breslow mice..")
  
  # Set cumulative hazards as missing (they will be updated)
  dat[, c("H1", "H2")] <- NA
  
  mat_breslow["X", c("Z", "eps", "H1", "H2")] <- 1
  meths_mice["H1"] <- paste("~I(", expression(update_basehaz(t, ev1, X, Z)),")")
  meths_mice["H2"] <- paste("~I(", expression(update_basehaz(t, ev2, X, Z)),")")
  
  imp_breslow <- mice::mice(
    dat, 
    m = m,
    method = meths_mice, 
    predictorMatrix = mat_breslow,
    maxit = 10, # should be enough
    print = FALSE
  )
  
  # With interaction
  dat[, c("H1_Z", "H2_Z")] <- NA
  mat_breslow["X", c("H1_Z", "H2_Z")] <- 1
  meths_mice["H1_Z"] <- "~I(H1 * Z)"
  meths_mice["H2_Z"] <- "~I(H2 * Z)"
  
  imp_breslow_int <- mice::mice(
    dat, 
    m = m,
    method = meths_mice, 
    predictorMatrix = mat_breslow,
    maxit = 10, # should be enough
    print = FALSE
  )
  
  # -- Summarise
  
  complist <- list(
    "imp_nels" = mice::complete(imp_nels, action = "all"),
    "imp_true" = mice::complete(imp_true, action = "all"),
    "imp_breslow" = mice::complete(imp_breslow, action = "all"),
    "imp_nels_int" = mice::complete(imp_nels_int, action = "all"),
    "imp_true_int" = mice::complete(imp_true_int, action = "all"),
    "imp_breslow_int" = mice::complete(imp_breslow_int, action = "all"),
    "smcfcs" = imp_smcfcs$value$impDatasets
  )
  
  mods_complist <- purrr::modify_depth(complist, .depth = 2, ~ setup_mstate(.x)) 
  
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
        analy == "imp_nels" ~ as.numeric(!is.null(imp_nels$loggedEvents)),
        analy == "imp_true" ~ as.numeric(!(is.null(imp_true$loggedEvents))),
        analy == "imp_breslow" ~ as.numeric(!(is.null(imp_breslow$loggedEvents))),
        analy == "imp_nels_int" ~ as.numeric(!is.null(imp_nels_int$loggedEvents)),
        analy == "imp_true_int" ~ as.numeric(!(is.null(imp_true_int$loggedEvents))),
        analy == "imp_breslow_int" ~ as.numeric(!(is.null(imp_breslow_int$loggedEvents))),
        analy %in% c("CCA", "ref") ~ 0
      ) 
    ) %>%
    cbind.data.frame(scenario, row.names = NULL)

  # Clear from memory
  rm(
    "imp_smcfcs", "imp_nels", "imp_true", "imp_breslow", 
    "imp_nels_int", "imp_true_int", "imp_breslow_int"
  )
  
  return(estimates)
}


# Run sims ----------------------------------------------------------------


# - save a list on shark, rbindlist locally (pull repo on shark)
# - add new functions to R/ in package


# Run 
scen_num <- Sys.getenv("SLURM_ARRAY_TASK_ID")
scenario <- scenarios_raw[scen_num, ]
n_sim <- scenario$n_sim

res <- lapply(
  X = seq_len(n_sim),
  FUN = function(rep) {
    one_simulation_breslow(
      scenario = scenario, 
      n_cpu = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")), 
      rep_num = rep
    )
  }
)

# Give a name depending on scenario..
res_filename <- paste0(
  "analysis/supplement-simulations/",
  "suppl-sims_breslow_scen", scen_num,
  ".rds"
)
saveRDS(res, file = res_filename)
