devtools::load_all()
library(survival)
library(mstate)
library(mice)
source("analysis/supplement-simulations/helpers-three-comprisks.R")
source("analysis/supplement-simulations/scenarios-supplemental.R")

one_simulation_threecomp <- function(scenario, n_cpu, rep_num) {
  
  # Set the seed
  scen_number <- as.numeric(rownames(scenario))
  set.seed(rep_num * scen_number)
  
  # Set up parameters
  baseline <- CauseSpecCovarMI::mds_shape_rates
  
  # Different basehaz setting, with another constant hazard as third
  ev1_pars <- list("a1" = 1.5, "h1_0" = 0.04, "b1" = scenario$beta1, "gamm1" = 1)
  
  ev2_pars <- list(
    "a2" = baseline[baseline$state == "NRM", "shape"], 
    "h2_0" = baseline[baseline$state == "NRM", "rate"], 
    "b2" = .5, 
    "gamm2" = .5
  )
  
  # constant; plot true baseline cumincs somewhere..
  ev3_pars <- list("a3" = 1, "h3_0" = 0.1, "b3" = 0.75, "gamm3" = .25)
  
  dat <- generate_dat_threecomp(
    n = scenario$n,
    X_type = scenario$X_level, 
    r = 0.5,
    ev1_pars, ev2_pars, ev3_pars,
    rate_cens = baseline[baseline$state == "EFS", "rate"],
    mech = scenario$miss_mech,
    eta1 = scenario$eta1,
    p = scenario$prop_miss
  )
  
  # -- Run full dataset and CCA
  
  mod_ref <- setup_mstate_threecomp(dat %>% dplyr::mutate(X = .data$X_orig))
  mod_CCA <- setup_mstate_threecomp(dat)
  
  # -- smcfcs
  
  m <- 50 # make 25 or 50
  meths_smcfcs <- mice::make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))
  
  imp_smcfcs <- record_warning(
    smcfcs::smcfcs.parallel(
      seed = round(runif(1, 0, 1e5)),# change somehow here; random based on initial seed
      n_cores = n_cpu,
      cl_type = "FORK",
      originaldata = dat, 
      smtype = "compet", 
      smformula = c(
        "Surv(t, eps == 1) ~ X + Z",
        "Surv(t, eps == 2) ~ X + Z",
        "Surv(t, eps == 3) ~ X + Z"
      ), 
      method = meths_smcfcs, 
      m = m, 
      numit = 20, # drop to 20 
      rjlimit = 5000
    ) # 5 times higher than default, avoid rej sampling errors
  )
  
  # -- mice setup
  
  mat <- mice::make.predictorMatrix(dat) 
  mat[] <- 0L 
  mat_ch123 <- mat_ch123_int <- mat 
  
  mat_ch123["X", c("Z", "eps", "H1", "H2", "H3")] <- 1
  mat_ch123_int["X", c("Z", "eps", "H1", "H2", "H3", "H1_Z", "H2_Z", "H3_Z")] <- 1
  meths_mice <-  mice::make.method(dat, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
  
  # -- mice standard and interaction
  
  # ch123
  imp_ch123 <- mice::mice(
    dat, 
    m = m,
    method = meths_mice, 
    predictorMatrix = mat_ch123,
    maxit = 1, 
    print = FALSE
  )
  
  # With interaction
  imp_ch123_int <- mice::mice(
    dat, 
    m = m,
    method = meths_mice, 
    predictorMatrix = mat_ch123_int,
    maxit = 1, 
    print = FALSE
  )
  
  # -- Summarise
  
  complist <- list(
    "imp_ch123" = mice::complete(imp_ch123, action = "all"),
    "imp_ch123_int" = mice::complete(imp_ch123_int, action = "all"),
    "smcfcs" = imp_smcfcs$value$impDatasets
  )
  
  mods_complist <- purrr::modify_depth(complist, .depth = 2, ~ setup_mstate_threecomp(.x)) 
  
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
        stringr::str_detect(var, "Z.2") ~ ev2_pars$gamm2,
        stringr::str_detect(var, "X.3") ~ ev3_pars$b3,
        stringr::str_detect(var, "Z.3") ~ ev3_pars$gamm3
      )
    ) %>% 
    
    # Add mice warnings and rej sampling errors for smcfcs
    dplyr::mutate(
      warns = dplyr::case_when(
        analy == "smcfcs" ~ as.numeric(stringr::str_extract(imp_smcfcs$warning, "[0-9]+")),
        analy == "imp_ch123" ~ as.numeric(!is.null(imp_ch123$loggedEvents)),
        analy == "imp_ch123_int" ~ as.numeric(!(is.null(imp_ch123_int$loggedEvents))),
        analy %in% c("CCA", "ref") ~ 0
      ) 
    ) %>%
    cbind.data.frame(scenario, row.names = NULL)
  
  # Clear from memory
  rm("imp_smcfcs", "imp_ch123", "imp_ch123_int")
  
  # Notify
  cat(paste0("\n", "End of rep #", as.character(rep_num), "!\n"))
  
  return(estimates)
}


# Run scenario ------------------------------------------------------------

# Run 
scen_num <- Sys.getenv("SLURM_ARRAY_TASK_ID")
scenario <- scenarios_raw[scen_num, ]
n_sim <- scenario$n_sim

res <- lapply(
  X = seq_len(n_sim),
  FUN = function(rep) {
    one_simulation_threecomp(
      scenario = scenario, 
      n_cpu = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")), 
      rep_num = rep
    )
  }
)

# Give a name depending on scenario..
res_filename <- paste0(
  "analysis/supplement-simulations/",
  "suppl-sims_threecomp_scen", scen_num,
  ".rds"
)
saveRDS(res, file = res_filename)
