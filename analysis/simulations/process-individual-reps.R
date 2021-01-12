

# Run this with large memory..

# Regression coefficients
sims_regr_full <- list.files(
  path = "data/sim-reps_indiv/regr/",
  full.names = TRUE
) %>% 
  pbapply::pblapply(readRDS) %>% 
  data.table::rbindlist()

# Save
# fst::write_fst(sim_regr_full, path = "data/sims_preds_full.fst")


# Predictions
sims_preds_full <- list.files(
  path = "data/sim-reps_indiv/preds/",
  full.names = TRUE
) %>% 
  pbapply::pblapply(readRDS) %>% 
  data.table::rbindlist()

# Save
# fst::write_fst(sims_preds_full, path = "data/sims_preds_full.fst")
