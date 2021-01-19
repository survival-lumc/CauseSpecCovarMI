##***************************##
## Generating synthetic data ##
##***************************##


# (For possible sharing..) - based on synthpop
options(contrasts = rep("contr.treatment", 2)) 

# Read-in (without admin censoring)
dat_mds <- fst::read_fst("data/dat-mds.fst") %>% 
  data.table::setDT()


<<<<<<< HEAD
# Prepare synthpop --------------------------------------------------------


dat_mds[, ':=' (
  srv_s_allo1 = NULL,
  srv_allo1 = NULL
)]

x <- dat_mds[, !c("ci_allo1", "ci_s_allo1")]
test_syn <- synthpop::syn(x)
xp <- test_syn$syn
synthpop::compare(test_syn, dat_mds)

naniar::gg_miss_upset(dat_mds)
naniar::gg_miss_upset(xp)

# Make comprisk syn function - based on riskreg? what about cases with missing covars?
syn.cmprskCSC <- function(y, x, xp, delta = dat_mds$ci_s_allo1, ...) {
  
  mod <- riskRegression::CSC(Hist(y, delta) ~ ., data = cbind.data.frame(y, delta, x))
  
  # Return
}

# Ideas: 
# 1. synthetise covars without outcome, keep missings
# 2. Impute m = 1  (enough iters) with smcfcs compatible with *original* outcome
# 3. Fit survreg on this new hybrid data, store baseline shape and rate (intercept) 
#   (We don't on original since we do not know how to pool shape?)
# 4. From imputations on original dataset, store the final PH model coefficients
# 5. Compute linear predictor using these coefficients, and hybrid dataset covariates
# 6. Generate latent event times
# 7. Add censoring also form weibull survreg


# Test new synth idea -----------------------------------------------------

outcome <- c("ci_allo1", "ci_s_allo1")

meths <- rep("cart", ncol(dat_mds))
names(meths) <- colnames(dat_mds)
meths[outcome] <- ""
meths[!(names(meths) %in% outcome)][1] <- "sample" 
meths

=======
# Prepare model formulas --------------------------------------------------


outcomes <- c("ci_s_allo1", "ci_allo1", "srv_s_allo1", "srv_allo1")
predictors <- sort(colnames(dat_mds)[!(colnames(dat_mds) %in% outcomes)])

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
  function(imp) survival::survreg(
    formula = Surv(ci_allo1, ci_s_allo1 == 0) ~ 1, 
    dist = "weibull", data = imp
  )
)  %>% 
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


outcomes <- c("ci_allo1", "ci_s_allo1")#

# We synth outcomes as cart, but they are use just for placeholder imp
meths <- rep("cart", ncol(dat_mds))
names(meths) <- colnames(dat_mds)
meths[c("srv_s_allo1", "srv_allo1")] <- ""

# Sort visiting sequence - from least missing to highest
miss_pct <- sapply(dat_mds, function(x) mean(is.na(x)))
visit_seq <- order(miss_pct)[which(
  !(names(meths) %in% c(outcomes, "srv_s_allo1", "srv_allo1"))
)]
meths[visit_seq[1]] <- "sample"
visit_seq <- c(visit_seq, which(names(meths) %in% outcomes))

# Check it
meths
visit_seq
>>>>>>> 81e2e71f2ca9e31cae01e31c8793c3ea8c400d90

synth_covar <- synthpop::syn(
  data = dat_mds,
  method = meths, 
<<<<<<< HEAD
  #visit.sequence = c(order(miss_props), which(meths == "")), 
  proper = TRUE, 
  models = TRUE,
  seed = 2021
)

# Use smfcs
predictors <- sort(
  colnames(synth_covar$syn)[!(colnames(synth_covar$syn) %in% outcome)]
) 

# Both REL and NRM have same rhs
rhs <- paste(predictors, collapse = " + ")

# Make both model formulas
form_rel <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 1) ~ ", rhs))
form_nrm <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", rhs))

var_names_miss <- colnames(synth_covar$syn)[sapply(synth_covar$syn, anyNA)]

meths <- set_mi_methods(
  dat = synth_covar$syn,
  var_names_miss = var_names_miss, 
=======
  visit.sequence = visit_seq,
  proper = TRUE, 
  seed = 2020,
  minnumlevels = 4 # to make ci_s_allo1 factor
)

# Compare
synthpop::compare(synth_covar, dat_mds)
naniar::gg_miss_upset(dat_mds)
naniar::gg_miss_upset(synth_covar$syn)

# Check any duplicates with original - non!
sum(duplicated(rbind(synth_covar$syn, dat_mds)))

# Impute syn covar once ---------------------------------------------------


meths <- set_mi_methods(
  dat = synth_covar$syn,
  var_names_miss = naniar::miss_var_which(synth_covar$syn), 
>>>>>>> 81e2e71f2ca9e31cae01e31c8793c3ea8c400d90
  imp_type = "smcfcs",
  cont_method = "norm"
)

# Make formula
smform_smcfcs <- c(
  Reduce(paste, deparse(form_rel)),
  Reduce(paste, deparse(form_nrm))
)

<<<<<<< HEAD
imps_smcfcs <- smcfcs::smcfcs(
  originaldata = synth_covar$syn,
  smtype = "compet",
  smformula = smform_smcfcs,
  m = 1,
  numit = 5, # make 25
  method = meths#,
  #rjlimit = 5000 # 5x normal
)
impo <- imps_smcfcs$impDatasets[[1]]

#plot(imps_smcfcs)
mod_rel <- survival::survreg(
  formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ .,
  data = impo
)
coefs_rel <- coefficients(mod_rel)
shape_rel <- 1/ mod_rel$scale
baserate_rel <- exp(coefs_rel[1])^(-shape_rel)
#instead use coefs from smcfcs
ph_coefs <- -coefs_rel[-1] * shape_rel 
mod_mat <- model.matrix(~., data = impo)[,-(1:3)]
lp_rel <- drop(mod_mat %*% ph_coefs) # here is the lp


# nrm ---------------------------------------------------------------------



mod_nrm <- survival::survreg(
  formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ .,
  data = impo
)
coefs_nrm <- coefficients(mod_nrm)
shape_nrm <- 1/ mod_nrm$scale
baserate_nrm <- exp(coefs_nrm[1])^(-shape_nrm)
#instead use coefs from smcfcs
ph_coefs_nrm <- -coefs_nrm[-1] * shape_nrm 
mod_mat <- model.matrix(~., data = impo)[,-(1:3)]
lp_nrm <- drop(mod_mat %*% ph_coefs_nrm) # here is the lp

# Need also censoring..

# Bind to syndata
newdato <- synth_covar$syn %>% 
  data.table() %>% 
  .[, ':=' (
    lp_rel = exp(lp_rel),
    lp_nrm = exp(lp_nrm),
    ci_allo1 = NULL,
    ci_s_allo1 = NULL
=======
# Impute 
imps_synth <- parlSMCFCS::parlsmcfcs(
  originaldata = synth_covar$syn,
  smtype = "compet",
  smformula = smform_smcfcs,
  m = 1, # should be 1
  numit = 20,
  seed = 2021,
  n_cores = 1,
  method = meths
)

# Remove synthesised outcomes
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
>>>>>>> 81e2e71f2ca9e31cae01e31c8793c3ea8c400d90
  )] %>% 
  
  # Gen new event times
  .[, ':=' (
<<<<<<< HEAD
    t_rel = rweibull_KM(.N, shape_rel, baserate_rel * lp_rel),
    t_nrm = rweibull_KM(.N, shape_nrm, baserate_nrm * lp_nrm)
  )] %>% 
  
  # 


# Rest --------------------------------------------------------------------





# How to factorise? Less to more missing
outcome <- c("ci_allo1", "ci_s_allo1")
miss_props <- unlist(lapply(dat_mds[, !..outcome], function(col) mean(is.na(col))))
sort(miss_props)

# Set methods
meths <- rep("cart", ncol(dat_mds))
names(meths) <- colnames(dat_mds)
meths[names(sort(miss_props)[1])] <- "sample" 

# Edit outcomes, event indicators synthesised with survtree
meths[c("ci_allo1")] <- "survctree" 
meths[c("ev1", "ev2")] <- "" 


# Create synthetic data ---------------------------------------------------


# Browser() will unfortunately open, just press continue
dat_mds_synth <- synthpop::syn(
  data = dat_mds,
  method = meths, 
  visit.sequence = c(order(miss_props), which(meths == "survctree")), 
  proper = TRUE, 
  models = TRUE,
  event = list("t_ev1" = "ev1", "t_ev2" = "ev2"),
  seed = 2021
)

# Compare
synthpop::compare(dat_mds_synth, dat_mds)

# Can compare model results later..
=======
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


# Compare cca -------------------------------------------------------------

# Read-in non admin cens one
dat_mds <- fst::read_fst("data/dat-mds_admin-cens.fst")

#
cbind(
  "real" = coef(survival::coxph(form_rel, data = dat_mds)),
  "synth" = coef(survival::coxph(form_rel, data = dat_mds_syn))
)

cbind(
  "real" = coef(survival::coxph(form_nrm, data = dat_mds)),
  "synth" = coef(survival::coxph(form_nrm, data = dat_mds_syn))
)
>>>>>>> 81e2e71f2ca9e31cae01e31c8793c3ea8c400d90


# Save synthetic data -----------------------------------------------------


<<<<<<< HEAD
fst::write_fst(x = dat_mds_synth$syn, path = "analysis/data/dat-mds-synth.fst")
=======
fst::write_fst(
  x = dat_mds_syn, 
  path = "data/dat-mds-synth.fst"
)
>>>>>>> 81e2e71f2ca9e31cae01e31c8793c3ea8c400d90
