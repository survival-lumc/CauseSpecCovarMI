##******************************##
## Summarise simulation results ##
##******************************##

# Process individual replications -----------------------------------------

# Regression coefficients
sims_regr_full <- list.files(
  path = "data/sim-reps_indiv/regr/",
  full.names = TRUE
) %>% 
  pbapply::pblapply(readRDS) %>% 
  data.table::rbindlist()

# Save
fst::write_fst(sim_regr_full, path = "data/sims_preds_full.fst")


# Predictions
sims_preds_full <- list.files(
  path = "data/sim-reps_indiv/preds/",
  full.names = TRUE
) %>% 
  pbapply::pblapply(readRDS) %>% 
  data.table::rbindlist()

# Save
fst::write_fst(sims_preds_full, path = "data/sims_preds_full.fst")

rm(sims_regr_full, sims_preds_full)


# Summarise regr ----------------------------------------------------------


# Read-in estimates
sims_regr_full <- fst::read_fst("data/sims_regr_full.fst") %>% 
  data.table::setDT() %>% 
  CauseSpecCovarMI::extract_scen_num()

# Compute measures per scenario
sims_regr_summary <- sims_regr_full[, .(
    n = .N,
    est = mean(estimate),
    se = mean(std.error),
    se_mcse = sqrt(var(std.error^2) / (4 * .N * mean(std.error^2))),
    emp_se = sd(estimate),
    cover = mean(`2.5 %` < true & true < `97.5 %`),
    bias = mean(estimate - true),
    rmse = sqrt(mean((estimate - true)^2)),
    rmse_mcse = CauseSpecCovarMI::rmse_mcse(estimate, true, .N),
    warns = mean(warns)
  ), by = .(var, m, analy, true, scen_summary, scen_num)] %>% 
  
  # Add mcse bias and coverage
  .[, ':=' (
    bias_mcse = emp_se / sqrt(n),
    cover_mcse = sqrt((cover * (1 - cover)) / n)
  )] %>% 
  CauseSpecCovarMI::format_scen_summary()

# Save sims
data.table::setorder(sims_regr_summary, var, analy, m)
fst::write_fst(x = sims_regr_summary, path = "data/sims_regr_summary.fst")

# Clear environment
rm(sims_regr_full)


# Summarise predictions ---------------------------------------------------


sims_preds_full <- fst::read_fst("data/sims_preds_full.fst") %>% 
  data.table::setDT() %>% 
  CauseSpecCovarMI::extract_scen_num()

# Compute measures per scenario
sims_preds_summary <- sims_preds_full[, .(
  n = .N,
  prob = mean(p_pool),
  emp_se = sd(p_pool),
  bias = mean(p_pool - true),
  rmse = sqrt(mean((p_pool - true)^2)),
  rmse_mcse = CauseSpecCovarMI::rmse_mcse(p_pool, true, .N)
), by = .(analy, m,`combo-X_Z`, times, state, true, scen_summary, scen_num)] %>% 
  
  # Add mcse bias and coverage
  .[, ':=' (
    bias_mcse = emp_se / sqrt(n)
  )] %>% 
  CauseSpecCovarMI::format_scen_summary() %>% 
  
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
data.table::setorder(sims_preds_summary, state, times, analy, m)
fst::write_fst(x = sims_preds_summary, path = "data/sims_preds_summary.fst")

#rm(sims_preds_full)
