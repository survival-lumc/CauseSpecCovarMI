##**********************##
## Analysis of mds data ##
##**********************##

library(CauseSpecCovarMI)

# Choose whether synthetic or not
synth <- FALSE

if (synth) {
  dat_mds <- CauseSpecCovarMI::dat_mds_synth %>% data.table::setDT()
} else {
  dat_mds <- fst::read_fst("data/dat-mds_admin-cens.fst") %>% data.table::setDT()
}

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2)) 


# Prepare model formulas --------------------------------------------------


outcomes <- c("ci_s_allo1", "ci_allo1", "srv_s_allo1", "srv_allo1")
predictors <- sort(colnames(dat_mds)[!(colnames(dat_mds) %in% outcomes)]) 

# Both REL and NRM have same rhs
rhs <- paste(predictors, collapse = " + ")

# Make both model formulas
form_rel <- stats::reformulate(
  termlabels = predictors,
  response = "Surv(ci_allo1, ci_s_allo1 == 1)"
)
  
form_nrm <- stats::reformulate(
  termlabels = predictors,
  response = "Surv(ci_allo1, ci_s_allo1 == 2)"
)

# Prepare MI --------------------------------------------------------------


# First, for competing risks, we make indicators for comp events and
# compute cumulative hazards
dat_mds[, ':=' (
  ev1 = ifelse(ci_s_allo1 == 1, 1, 0),
  ev2 = ifelse(ci_s_allo1 == 2, 1, 0)
)]

dat_mds[, ':=' (
  H1 = CauseSpecCovarMI::nelsaalen_timefixed(dat = data.frame(.SD), timevar = ci_allo1, statusvar = ev1),
  H2 = CauseSpecCovarMI::nelsaalen_timefixed(dat = data.frame(.SD), timevar = ci_allo1, statusvar = ev2)
)]

# Which vars actually have some data missing? (no binary vars)
var_names_miss <- colnames(dat_mds)[sapply(dat_mds, anyNA)]

# Set methods accordingly
meths <- CauseSpecCovarMI:::set_mi_methods(
  dat = dat_mds,
  var_names_miss = var_names_miss, 
  imp_type = "mice",
  cont_method = "norm"
)

# Create matrix for MI
matpred <- matrix(
  data = 1, 
  ncol = ncol(dat_mds), 
  nrow = ncol(dat_mds),
  dimnames = list(names(dat_mds), names(dat_mds))
)

# Don't impute a var using itself + use only H1 and H2 as outcomes
# + don't impute any complete variable
diag(matpred) <- 0 
matpred[, outcomes] <- 0
matpred[!(rownames(matpred) %in% var_names_miss), ] <- 0
matpred

# Set number of imputations and iterations globally
m <- 2 #100
iters <- 2 #20 


# Run MICE ----------------------------------------------------------------


# Try in parallel
imps_mice <- mice::parlmice(
  data = dat_mds,
  m = m,
  maxit = iters,
  cluster.seed = 1984, 
  n.core = 2, #10, 
  n.imp.core = 1, #10, 
  method = meths,
  predictorMatrix = matpred
)

# save imputations..
# saveRDS(imps_mice, file = "data/imps-mds-mice.rds")
# rm(imps_mice)


# Run smcfcs --------------------------------------------------------------


# Change coding of polr and polyreg, exclude the prev vars
meths <- CauseSpecCovarMI::set_mi_methods(
  dat = dat_mds,
  var_names_miss = var_names_miss, 
  imp_type = "smcfcs",
  cont_method = "norm"
)

# Make formula
smform_smcfcs <- c(
  Reduce(paste, deparse(form_rel)),
  Reduce(paste, deparse(form_nrm))
)

imps_smcfcs <- parlSMCFCS::parlsmcfcs(
  seed = 4891,
  n_cores = 2, #10,
  outfile = "out_paralsmcfcs.txt",
  originaldata = dat_mds,
  smtype = "compet",
  smformula = smform_smcfcs,
  m = 2, #m,
  numit = 2, #iters,
  method = meths,
  rjlimit = 3000 # 3x normal, reduce num warnings
)


# save imputations as .rds..
#saveRDS(imps_smcfcs, file = "data/imps-mds-smcfcs.rds")
#rm(imps_smcfcs)


# Pool regression coefficients --------------------------------------------


# Prepare lists - read-in imputations
imps_mice <- readRDS("data/imps-mds-mice.rds")
imps_smcfcs <- readRDS("data/imps-mds-smcfcs.rds")
impdats_mice <- mice::complete(imps_mice, action = "all")
impdats_smcfcs <- imps_smcfcs$impDatasets


# Can verify n imps with howManyImputations::how_many_imputations()
# on object left after pool()

# Relapse models ----------------------------------------------------------

mice_rel <- pbapply::pblapply(
  impdats_mice, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE, exponentiate = TRUE)

smcfcs_rel <- pbapply::pblapply(
  impdats_smcfcs, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE, exponentiate = TRUE)


# NRM models --------------------------------------------------------------


mice_nrm <- pbapply::pblapply(
  impdats_mice, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE, exponentiate = TRUE)

# smcfcs
smcfcs_nrm <- pbapply::pblapply(
  impdats_smcfcs, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE, exponentiate = TRUE)


# CCA ---------------------------------------------------------------------

  
mod_rel <- survival::coxph(form_rel, data = dat_mds) %>% 
  summary(conf.int = 0.95) %$% 
  conf.int %>% 
  data.table::data.table(keep.rownames = "term") %>% 
  data.table::setnames(
    old = c("exp(coef)", "lower .95", "upper .95"), 
    new = c("estimate", "2.5 %", "97.5 %")
  )

mod_nrm <- survival::coxph(form_nrm, data = dat_mds) %>% 
  summary(conf.int = 0.95) %$% 
  conf.int %>% 
  data.table::data.table(keep.rownames = "term") %>% 
  data.table::setnames(
    old = c("exp(coef)", "lower .95", "upper .95"), 
    new = c("estimate", "2.5 %", "97.5 %")
  )

# Check proportionality of complete cases
survival::cox.zph(survival::coxph(form_rel, data = dat_mds))
survival::cox.zph(survival::coxph(form_nrm, data = dat_mds))

# Look into one of imputed datasets
# - crnocr and patient age test as non-proportional
# - residuals show it is likely not too dramatic
# Note: better method for checking proportionality with MI is described
# in Keogh et al. (2018), using joint wald tests
zph_rel <- survival::cox.zph(survival::coxph(form_rel, data = impdats_mice[[15]]))
zph_nrm <- survival::cox.zph(survival::coxph(form_nrm, data = impdats_mice[[15]]))
plot(zph_rel)
plot(zph_nrm)


# Forest plot -------------------------------------------------------------


# Read in data dictionary
dictionary <- get(load("R/data_dictionary.rda"))

list_rel <- list(
  "SMC-FCS" = smcfcs_rel,
  "MICE" = mice_rel,
  "CCA" = mod_rel
)

list_nrm <- list(
  "SMC-FCS" = smcfcs_nrm,
  "MICE" = mice_nrm,
  "CCA" = mod_nrm
)

forest_rel <- CauseSpecCovarMI:::ggplot_grouped_forest(
  dat = dat_mds,
  dictionary = dictionary,
  results = list_rel,
  event = "Relapse",
  form = form_rel
)

forest_nrm <- CauseSpecCovarMI:::ggplot_grouped_forest(
  dat = dat_mds,
  dictionary = dictionary,
  results = list_nrm,
  event = "Non-relapse mortality",
  form = form_nrm,
  lims_x = c(0.5, 3.1)
)

# Change to rel
ggplot2::ggsave(
  filename = "analysis/figures/forest_rel.eps", 
  plot = forest_rel, 
  width = 10, 
  height = 11
)

# Change to nrm
ggplot2::ggsave(
  filename = "analysis/figures/forest_nrm.eps", 
  plot = forest_rel, 
  width = 10, 
  height = 11
)


# Predictions -------------------------------------------------------------


# Prepare tmat and msprep
tmat <- mstate::trans.comprisk(2, c("Rel", "NRM"))
dat_mds[, ':=' (ev1 = ci_s_allo1 == 1, ev2 = ci_s_allo1 == 2)]

dat_msp_mds <- mstate::msprep(
  time = c(NA, "ci_allo1", "ci_allo1"),
  status = c(NA, "ev1", "ev2"), 
  data = data.frame(dat_mds),
  trans = tmat,
  keep = predictors
) 

# Expand covariates 
dat_expand_mds <- mstate::expand.covs(
  dat_msp_mds, predictors, append = TRUE, longnames = FALSE
)

# RHS formula mstate - keep transition-specific coefs
form_mds_mstate <- reformulate(
  termlabels = c(
    colnames(dat_expand_mds)[grep(x = colnames(dat_expand_mds), pattern = "\\.[1-9]+$")],
    "strata(trans)"
  ),
  response = "Surv(time, status)"
)

rm(dat_expand_mds)

# Choose patients to predict
ref_pat <- dat_mds[, lapply(.SD, function(col) {
  CauseSpecCovarMI::choose_standard_refpat(col, "median", "reference") 
}), .SD = predictors]

# Make other ref_pats - for excess blasts and sAML
ref_pat_excess <- data.table::copy(ref_pat)
ref_pat_excess[1, mdsclass := "MDS with excess blasts"]

ref_pat_saml <- data.table::copy(ref_pat)
ref_pat_saml[1, mdsclass := "sAML"]

#newpat <- make_mstate_refpat(ref_pat, tmat, predictors) 

ref_pats <- list(
  "MDS without excess blasts" = CauseSpecCovarMI::make_mstate_refpat(ref_pat, tmat, predictors),
  "MDS with excess blasts" = CauseSpecCovarMI::make_mstate_refpat(ref_pat_excess, tmat, predictors),
  "sAML" = CauseSpecCovarMI::make_mstate_refpat(ref_pat_saml, tmat, predictors)
)

rm(ref_pat, ref_pat_excess, ref_pat_saml)


# Predict for mice and smcfcs - do for all later
preds_list_mice <- pbapply::pblapply(impdats_mice[1:5], function(impdat) {
  
  # Run model
  mod <- CauseSpecCovarMI:::run_mds_model(
    form = form_mds_mstate,
    tmat = tmat,
    dat = impdat
  )
  
  # Predict
  preds <- CauseSpecCovarMI:::predict_mds_model(
    mod = mod, 
    ref_pats = ref_pats,
    horizon = 5, 
    tmat = tmat 
  )
  
  return(preds)
})

preds_list_smcfcs <- pbapply::pblapply(impdats_smcfcs[1:5], function(impdat) {
  
  # Run model
  mod <- CauseSpecCovarMI:::run_mds_model(
    form = form_mds_mstate,
    tmat = tmat,
    dat = impdat
  )
  
  # Predict
  preds <- CauseSpecCovarMI:::predict_mds_model(
    mod, 
    ref_pats = ref_pats,
    horizon = 5, 
    tmat = tmat 
  )
  
  return(preds)
})

# Pool both
preds_mice <- preds_list_mice %>% CauseSpecCovarMI::pool_morisot(by_vars = c("state", "ref_pat"))
preds_smcfcs <- preds_list_smcfcs %>% CauseSpecCovarMI::pool_morisot(by_vars = c("state", "ref_pat"))

# Make predictions for CCA
preds_CCA <- CauseSpecCovarMI::predict_mds_model(
  mod = run_mds_model(form_mds_mstate, tmat, dat_mds), 
  ref_pats = ref_pats,
  horizon = 5, 
  tmat = tmat
) %>% data.table::setDT()

# Recall transformation functions
cloglog <- Vectorize(function(x) log(-log(1 - x)))
inv_cloglog <- Vectorize(function(x) 1 - exp(-exp(x)))

# Add CI on invcloglog also for CCA
preds_CCA[, var_delta := se^2 / (log(1 - prob) * (1 - prob))^2]
preds_CCA[, ':=' (
  CI_low = inv_cloglog(cloglog(prob) - qnorm(0.975) * sqrt(var_delta)),
  CI_upp = inv_cloglog(cloglog(prob) + qnorm(0.975) * sqrt(var_delta)),
  p_pooled = prob # Just to bind rows later
)]

# Bind results together
preds <- rbind(preds_CCA, preds_mice, preds_smcfcs, fill = TRUE, idcol = "method")
preds[, method := factor(method, levels = 1:3, labels = c("$CC$", "$CH_{12}$", "smcfcs"))]
preds[, est := paste0(
  round(p_pooled * 100, 1), " [",
  round(CI_low * 100, 1), "; ",
  round(CI_upp * 100, 1), "]"
)]

# Present results for REL and NRM
preds_table <- data.table::dcast(
  preds[state != 1], 
  formula = state + ref_pat ~ method, value.var = "est"
)

preds_table[, ':=' (
  "MDS class" = factor(ref_pat, levels = c("MDS without excess blasts", "MDS with excess blasts", "sAML")),
  state = factor(state, levels = c(2, 3), labels = c("REL", "NRM"))
)]
data.table::setorder(preds_table, state, `MDS class`)

caption <- "Predicted cumulative incidence (\\%) at 5-years for reference patients with different MDS classes. 95\\% confidence intervals were constructed based on a complementary log-log transformation."

# library(kableExtra) add to suggests
kableExtra::kbl(
  x = preds_table[, c("MDS class", "$CC$", "$CH_{12}$", "smcfcs")], 
  format = "latex",
  booktabs = T,
  position = "ht",
  caption = caption,
  linesep = "",
  escape = F,
  digits = 1
) %>% 
  kableExtra::pack_rows(index = c("REL" = 3, "NRM" = 3))
