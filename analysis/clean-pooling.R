# Read imps
imps_smcfcs <- readRDS("data/imps-mds-smcfcs.rds")
impdats_smcfcs <- imps_smcfcs$impDatasets

# Read-in origin and prepare mstate
dat_mds <- fst::read_fst("data/dat-mds_admin-cens.fst") %>% 
  data.table::setDT()


# Prepare mstate ----------------------------------------------------------


outcomes_mds <- c("ci_s_allo1", "ci_allo1", "srv_s_allo1", "srv_allo1")
predictors_mds <- sort(colnames(dat_mds)[!(colnames(dat_mds) %in% outcomes_mds)]) 

# Prepare tmat and msprep
tmat <- mstate::trans.comprisk(2, c("Rel", "NRM"))
dat_mds[, ':=' (ev1 = ci_s_allo1 == 1, ev2 = ci_s_allo1 == 2)]

dat_msp_mds <- mstate::msprep(
  time = c(NA, "ci_allo1", "ci_allo1"),
  status = c(NA, "ev1", "ev2"), 
  data = as.data.frame(dat_mds),
  trans = tmat,
  keep = predictors_mds
) 

# Expand covariates 
dat_expand_mds <- mstate::expand.covs(
  dat_msp_mds, predictors_mds, append = TRUE, longnames = FALSE
)

# RHS formula mstate - keep transition-specific coefs
rhs_mds_mstate <- paste(
  colnames(dat_expand_mds)[grep(x = colnames(dat_expand_mds), pattern = "\\.[1-9]+$")],
  collapse = " + "
)

# Final formula
form_mds_mstate <- as.formula(
  paste0("Surv(time, status) ~ ", rhs_mds_mstate, " + strata(trans)")
)

# Remove redundant objects
rm(rhs_mds_mstate, dat_msp_mds)



# Refpat ------------------------------------------------------------------



# Make reference pat here 
ref_pat <- dat_mds[, lapply(.SD, function(col) {
  choose_standard_refpat(col, "median", "reference") 
}), .SD = predictors_mds]

# Make other ref_pat
ref_pat2 <- data.table::copy(ref_pat)
ref_pat2[1, mdsclass := "MDS with excess blasts"]

newpat <- make_mstate_refpat(ref_pat, tmat, predictors_mds) 

refos <- list(
  "MDS without excess blasts" = make_mstate_refpat(ref_pat, tmat, predictors_mds),
  "MDS with excess blasts" = make_mstate_refpat(ref_pat2, tmat, predictors_mds) #,
  #"sAML" = 
)

refos




# Make preds lists, in parallel?
# cl <- parallel::makeCluster(3, type = "PSOCK")
# parallel::clusterExport(cl, varlist = c("run_mds_model", "predict_mds_model",
#                                     "tmat", "refos", "form_mds_mstate"),
#                         envir = environment())
# parallel::clusterEvalQ(cl, expr = {
#   library(mstate)
#   library(data.table)
# })

preds_list <- pbapply::pblapply(impdats_smcfcs, function(impdat) {
  
  # Run model
  mod <- run_mds_model(
    form = form_mds_mstate,
    tmat = tmat,
    dat = impdat
  )
  
  # Predict
  preds <- predict_mds_model(
    mod, 
    ref_pats = refos,
    horizon = 5, 
    tmat = tmat 
  )
  
  return(preds)
}, cl = cl)

parallel::stopCluster(cl)

preds_list %>% 
  pool_morisot(by_vars = c("state", "ref_pat"))



