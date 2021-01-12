##*******************************##
## Summarise raw simulation data ##
##*******************************##


# Add these functions to .R
extract_scen_num <- function(dat) {
  dat[, ':=' (
    scen_num = gsub(
      pattern = ".*(scen_num=)|(-rep).*$", 
      replacement = "", 
      x = scen_summary
    ),
    scen_summary = gsub(
      pattern = "(-scen_num).*$", 
      replacement = "", 
      x = scen_summary
    )
  )]
  
  return(dat)
}

format_scen_summary <- function(dat) {
  
  # Labels for variables in the "scen_summary" column
  labs_scens <- c(
    "n", 
    "prop_miss", 
    "beta1", 
    "miss_mech", 
    "X_level",
    "rho", 
    "eta1", 
    "haz_shape"
  )
  
  # Adjust for eta minus label
  dat[, scen_summary := gsub(
    pattern = "eta1=-", 
    replacement = "eta1=min", 
    x = scen_summary
  )]
    
  # Separate the scen_summary column
  dat[, (labs_scens) := data.table::tstrsplit(scen_summary, split = "-")] 
    
  # Relevel new variables - ready for analysis
  dat[, ':=' (
    analy = factor(
      analy, levels = c("ref", "CCA", "ch1", "ch12", "ch12_int", "smcfcs")
    ),
    prop_miss = factor(
      prop_miss, 
      levels = c("prop_miss=0.1","prop_miss=0.5"),
      labels = c("10%", "50%")
    ),
    haz_shape = factor(
      haz_shape,
      levels = c("haz_shape=similar", "haz_shape=different"),
      labels = c("similar", "different")
    ),
    beta1 = factor(
      beta1,
      levels = c("beta1=0", "beta1=0.5", "beta1=1"),
      labels = c("0", "0.5", "1")
    ),
    eta1 = factor(
      eta1,
      levels = c("eta1=NA", "eta1=min1", "eta1=min2"),
      labels = c("None", "Weak", "Strong")
    ),
    miss_mech = factor(
      miss_mech,
      levels = c("miss_mech=MCAR", "miss_mech=MAR", 
                 "miss_mech=MAR_GEN", "miss_mech=MNAR"),
      labels = c("MCAR", "MAR", "MAR_GEN", "MNAR") # change to MAR-T ?
    ),
    X_level = factor(
      X_level,
      levels = c("X_level=continous", "X_level=binary"),
      labels = c("continuous", "binary")
    ),
    rho = factor(
      rho,
      levels = "rho=0.5",
      labels = "0.5"
    ),
    n = factor(
      n,
      levels = c("n=500", "n=2000"),
      labels = c("500", "2000")
    ),
    
    # Remove original summary variable
    scen_summary = NULL
  )]
  
  return(dat)
}

# Adapted from calc_absolute - without dplyr::pull for speed,
# operating solely on vectors
# https://github.com/meghapsimatrix/simhelpers/blob/master/R/calc_absolute.R
# https://cran.r-project.org/web/packages/simhelpers/vignettes/MCSE.html
rmse_mcse <- function(estimates, true, K) {
  
  # Keep first true value
  true_param <- true[1]

  # Calculate elements of rmse mcse
  t_bar <- mean(estimates) 
  var_t <- var(estimates) 
  t_bar_j <- (1 / (K - 1)) * (K * t_bar - estimates) 
  bias_j_sq <- (t_bar_j - true_param)^2 
  s_sq_t_j <- (1 / (K - 2)) * ((K - 1) * var_t - (K / (K - 1)) * (estimates - t_bar)^2) 
  rmse_j <- sqrt(bias_j_sq + s_sq_t_j) 
  mse <- mean((estimates - true_param)^2) 
  
  # Calculate rmse and mcse
  rmse <- sqrt(mse)
  mcse <- sqrt((1/(K)) * sum((rmse_j - rmse)^2))
  
  return(mcse)
}



# Regression coefficients -------------------------------------------------



# Read-in estimates
sims_regr_full <- fst::read_fst("data/sims_regr_full.fst") %>% 
  setDT() %>% 
  extract_scen_num()

# Compute measures per scenario
sims_regr_summary <- sims_regr_full[, .(
    n = .N,
    est = mean(estimate),
    se = mean(std.error),
    emp_se = sd(estimate),
    cover = mean(`2.5 %` < true & true < `97.5 %`),
    bias = mean(estimate - true),
    rmse = sqrt(mean((estimate - true)^2)),
    rmse_mcse = rmse_mcse(estimate, true, .N),
    warns = mean(warns)
  ), by = .(var, m, analy, true, scen_summary, scen_num)] %>% 
  
  # Add mcse bias and coverage
  .[, ':=' (
    bias_mcse = emp_se / sqrt(n),
    cover_mcse = sqrt((cover * (1 - cover)) / n)
  )] %>% 
  format_scen_summary()

# Save sims
setorder(sims_regr_summary, var, analy, m)
fst::write_fst(x = sims_regr_summary, path = "data/sims_regr_summary.fst")

# Clear environment
rm(sims_regr_full, sims_regr_summary)


# Predictions -------------------------------------------------------------

sims_preds_full <- fst::read_fst("data/sims_preds_full.fst") %>% 
  setDT() %>% 
  extract_scen_num()

# Compute measures per scenario
sims_preds_summary <- sims_preds_full[, .(
  n = .N,
  prob = mean(p_pool),
  emp_se = sd(p_pool),
  bias = mean(p_pool - true),
  rmse = sqrt(mean((p_pool - true)^2)),
  rmse_mcse = rmse_mcse(p_pool, true, .N)
), by = .(analy, m,`combo-X_Z`, times, state, true, scen_summary, scen_num)] %>% 
  
  # Add mcse bias and coverage
  .[, ':=' (
    bias_mcse = emp_se / sqrt(n)
  )] %>% 
  format_scen_summary() %>% 
  
  # Relevel prediction-specific factors
  .[, ':=' (
    state = factor(
      state,
      levels = c("1", "2", "3"),
      labels = c("EFS", "REL", "NRM")
    ),
    times = factor(
      times, 
      levels = c(0.5, 5, 10),
      labels = c("6 months", "5 years", "10 years")
    )
  )]

# Save sims
setorder(sims_preds_summary, state, times, analy, m)
fst::write_fst(x = sims_preds_summary, path = "data/sims_preds_summary.fst")

rm(sims_preds_full, sims_preds_summary)

