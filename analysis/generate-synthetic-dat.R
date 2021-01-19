##***************************##
## Generating synthetic data ##
##***************************##


# (For possible sharing..) - based on synthpop
options(contrasts = rep("contr.treatment", 2)) 

# Read-in (without admin censoring)
dat_mds <- fst::read_fst("data/dat-mds.fst") %>% 
  data.table::setDT()


# Prepare model formulas --------------------------------------------------


outcomes <- c("ci_s_allo1", "ci_allo1", "srv_s_allo1", "srv_allo1")
predictors <- sort(colnames(dat_mds)[!(colnames(dat_mds) %in% outcomes)]) # change order maybe later

# Both REL and NRM have same rhs
rhs <- paste(predictors, collapse = " + ")

# Make both model formulas
form_rel <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 1) ~ ", rhs))
form_nrm <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", rhs))



# Estimate Weibull parameters on original (imputed) data ------------------


# We use the first imputed dataset from smcfcs
imps_smcfcs <- readRDS("data/imps-mds-smcfcs.rds")

# Estimate relapse parameters
aft_rel_smcfcs <- lapply(
  imps_smcfcs$impDatasets, 
  function(imp) survival::survreg(
    formula = form_rel, dist = "weibull", data = imp)
) %>% 
  mice::pool() %>% 
  summary()


scale_rel <- exp(aft_rel_smcfcs[aft_rel_smcfcs$term == "Log(scale)", "estimate"])
intercept_rel <- aft_rel_smcfcs[aft_rel_smcfcs$term == "(Intercept)", "estimate"]
shape_rel <- 1/ scale_rel
baserate_rel <- exp(intercept_rel)^(-shape_rel)


# Estimate NRM parameters
aft_nrm_smcfcs <- lapply(
  imps_smcfcs$impDatasets, 
  function(imp) survival::survreg(
    formula = form_nrm, dist = "weibull", data = imp)
) %>% 
  mice::pool() %>% 
  summary()


scale_nrm <- exp(aft_nrm_smcfcs[aft_nrm_smcfcs$term == "Log(scale)", "estimate"])
intercept_nrm <- aft_nrm_smcfcs[aft_nrm_smcfcs$term == "(Intercept)", "estimate"]
shape_nrm <- 1/ scale_nrm
baserate_nrm <- exp(intercept_nrm)^(-shape_nrm)



#  Censoring --------------------------------------------------------------




# Estimate NRM parameters
aft_efs_smcfcs <- lapply(
  imps_smcfcs$impDatasets, 
  function(imp) {
    
    # Deal first with admin censoring
    imp[ci_allo1 == 10, ':=' (
      ci_allo1 = 10.1,
      ci_s_allo1 = 1 # Event just after
    )]
    
    # Run survreg model
    survival::survreg(
      formula = Surv(ci_allo1, ci_s_allo1 == 0) ~ 1, dist = "weibull", data = imp
    )
  })  %>% 
  mice::pool() %>% 
  summary()


scale_efs <- exp(aft_efs_smcfcs[aft_efs_smcfcs$term == "Log(scale)", "estimate"])
intercept_efs <- aft_efs_smcfcs[aft_efs_smcfcs$term == "(Intercept)", "estimate"]
shape_efs <- 1/ scale_efs
baserate_efs <- exp(intercept_efs)^(-shape_efs)


# Get pooled Cox coefficients for smcfcs and NRM  -------------------------


# These are the ones reported in the manuscript...
smcfcs_rel <- lapply(
  imps_smcfcs$impDatasets, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary()

smcfcs_nrm <- lapply(
  imps_smcfcs$impDatasets, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary()


# Prepare synthpop --------------------------------------------------------


outcomes <- c("ci_allo1", "ci_s_allo1", "srv_allo1", "srv_s_allo1")#

# Define methods - outcomes are not synthesised
meths <- rep("cart", ncol(dat_mds))
names(meths) <- colnames(dat_mds)
meths[outcomes] <- ""

# Sort visiting sequence
miss_pct <- sapply(dat_mds, function(x) mean(is.na(x)))
visit_seq <- order(miss_pct)[which(!(names(meths) %in% outcomes))]
meths[visit_seq[1]] <- "sample"
meths

synth_covar <- synthpop::syn(
  data = dat_mds,
  method = meths, 
  visit.sequence = visit_seq,
  proper = TRUE, 
  seed = 2021
)


# Compare
synthpop::compare(synth_covar, dat_mds)
naniar::gg_miss_upset(dat_mds)
naniar::gg_miss_upset(synth_covar$syn)



# Impute syn covar once ---------------------------------------------------


meths <- set_mi_methods(
  dat = synth_covar$syn,
  var_names_miss = naniar::miss_var_which(synth_covar$syn), 
  imp_type = "smcfcs",
  cont_method = "norm"
)

# Make formula
smform_smcfcs <- c(
  Reduce(paste, deparse(form_rel)),
  Reduce(paste, deparse(form_nrm))
)

# Impute 
imps_synth <- smcfcs::smcfcs(
  originaldata = synth_covar$syn,
  smtype = "compet",
  smformula = smform_smcfcs,
  m = 1,
  numit = 5, # make 25
  method = meths
)

# Remove outcomes
impdat_synth <- imps_synth$impDatasets[[1]] %>% 
  data.table::data.table() %>% 
  .[, !..outcomes]

# Make model matrix and compute linear predictors
mod_mat <- model.matrix(
  as.formula(paste0("~ ", rhs)), 
  data = impdat_synth
)[, -1]

rate_rel <- baserate_rel * exp(drop(mod_mat %*% smcfcs_rel$estimate))
rate_nrm <- baserate_nrm * exp(drop(mod_mat %*% smcfcs_nrm$estimate))

# Draw!
impdat_synth %>% 
  .[, ':=' (
    rate_rel = rate_rel,
    rate_nrm = rate_nrm
  )] %>% 
  
  # Gen new event times
  .[, ':=' (
    t_rel = mapply(rweibull_KM, rate_rel, MoreArgs = list(n = 1, alph = shape_rel)),
    t_nrm = mapply(rweibull_KM, rate_nrm, MoreArgs = list(n = 1, alph = shape_nrm)),
    t_cens = rweibull_KM(.N, shape_efs, baserate_efs) #
  )] %>% 
  
  # Make latent cmprsk time 
  .[, ':=' (
    t_tilde = pmin(t_rel, t_nrm),
    eps_tilde = ifelse(t_rel < t_nrm, 1, 2)
  )]  %>% 
  
  # Make final outcome
  .[, ':=' (
    ci_allo1 = pmin(t_tilde, t_cens),
    ci_s_allo1 = ifelse(t_cens < t_tilde, 0, eps_tilde)
  )] %>%

  # Admin censoring
  .[, ':=' (
    ci_allo1 = ifelse(ci_allo1 > 10, 10, ci_allo1),
    ci_s_allo1 = ifelse(ci_allo1 > 10, 0, ci_s_allo1)
  )]


# Bind these outcome to the synthesised one with missing 
dat_mds_syn <- synth_covar$syn %>% 
  data.table::data.table() %>% 
  .[, ':=' (
    ci_s_allo1 = impdat_synth$ci_s_allo1,
    ci_allo1 = impdat_synth$ci_allo1
  )]

# Try imputing
meths <- set_mi_methods(
  dat = dat_mds_syn,
  var_names_miss = naniar::miss_var_which(dat_mds_syn), 
  imp_type = "smcfcs",
  cont_method = "norm"
)

# Make formula
smform_smcfcs <- c(
  Reduce(paste, deparse(form_rel)),
  Reduce(paste, deparse(form_nrm))
)

# Impute 
imps_syn_full <- parlSMCFCS::parlsmcfcs(
  originaldata = dat_mds_syn,
  smtype = "compet",
  smformula = smform_smcfcs,
  m = 5,
  n_cores = 3,
  method = meths
)


syn_rel <- lapply(
  imps_syn_full$impDatasets, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary()

syn_rel$estimate
smcfcs_rel$estimate

# Compare cca -------------------------------------------------------------


table(impdat_synth$ci_s_allo1)
dat_mds <- fst::read_fst("data/dat-mds_admin-cens.fst")
table(dat_mds$ci_s_allo1)

survival::coxph(form_rel, data = dat_mds)
survival::coxph(form_rel, data = dat_mds_syn)


# Ideas: 
# 1. synthetise covars without outcome, keep missings
# 2. Impute m = 1  (enough iters) with smcfcs compatible with *original* outcome
# 3. Fit survreg on this new hybrid data, store baseline shape and rate (intercept) 
#   (We don't on original since we do not know how to pool shape?)
# 4. From imputations on original dataset, store the final PH model coefficients
# 5. Compute linear predictor using these coefficients, and hybrid dataset covariates
# 6. Generate latent event times
# 7. Add censoring also form weibull survreg


 


# Save synthetic data -----------------------------------------------------


fst::write_fst(x = dat_mds_synth$syn, path = "analysis/data/dat-mds-synth.fst")
