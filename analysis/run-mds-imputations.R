##**********************##
## Analysis of mds data ##
##**********************##


# Read-in
dat_mds <- fst::read_fst("data/dat-mds_admin-cens.fst") %>% 
  data.table::setDT()

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2)) 


# Prepare model formulas --------------------------------------------------


outcomes <- c("ci_s_allo1", "ci_allo1", "srv_s_allo1", "srv_allo1")
predictors <- sort(colnames(dat_mds)[!(colnames(dat_mds) %in% outcomes)]) # change order maybe later

# Both REL and NRM have same rhs
rhs <- paste(predictors, collapse = " + ")

# Make both model formulas
form_rel <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 1) ~ ", rhs))
form_nrm <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", rhs))


# Complete case analysis --------------------------------------------------


# Use mstate instead? or just for predictions??
cca_rel <- survival::coxph(form_rel, data = dat_mds)
cca_nrm <- survival::coxph(form_nrm, data = dat_mds)


# Prepare MI --------------------------------------------------------------


# First, for competing risks, we make indicators for comp events and
# compute cumulative hazards
dat_mds[, ':=' (
  ev1 = ifelse(ci_s_allo1 == 1, 1, 0),
  ev2 = ifelse(ci_s_allo1 == 2, 1, 0)
)]

dat_mds[, ':=' (
  H1 = nelsaalen_timefixed(dat = data.frame(.SD), timevar = ci_allo1, statusvar = ev1),
  H2 = nelsaalen_timefixed(dat = data.frame(.SD), timevar = ci_allo1, statusvar = ev2)
)]

# Which vars actually have some data missing? (no binary vars)
var_names_miss <- colnames(dat_mds)[sapply(dat_mds, anyNA)]

# Set methods accordingly
meths <- set_mi_methods(
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
m <- 5
iters <- 2 # test with 25


# Run MICE ----------------------------------------------------------------


# Try in parallel
imps_mice <- mice::parlmice(
  data = dat_mds,
  m = m,
  maxit = iters,
  cluster.seed = 1984, 
  n.core = 4, 
  n.imp.core = 1, #floor(m / 3),
  method = meths,
  predictorMatrix = matpred
)

# save imputations..
saveRDS(imps_mice, file = "data/imps-mds-mice.rds")


# Run smcfcs --------------------------------------------------------------


# Change coding of polr and polyreg, exclude the prev vars
meths <- set_mi_methods(
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
  n_cores = 4,
  outfile = "out_paralsmcfcs.txt",
  originaldata = dat_mds,
  smtype = "compet",
  smformula = smform_smcfcs,
  m = m,
  numit = iters,
  method = meths
  #rjlimit = 2000 # 2x normal, for donor age
)


# save imputations as .rds..
saveRDS(imps_smcfcs, file = "data/imps-mds-smcfcs.rds")


# Pool regression coefficients --------------------------------------------


# Prepare lists 
impdats_mice <- mice::complete(imps_mice, action = "all")
impdats_smcfcs <- imps_smcfcs$impDatasets

# Mice models
mice_rel <- lapply(
  impdats_mice, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)

howManyImputations::how_many_imputations(lapply(
  impdats_mice, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool())

mice_nrm <- lapply(
  impdats_mice, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)


# smcfcs models
smcfcs_rel <- lapply(
  impdats_smcfcs, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)

smcfcs_nrm <- lapply(
  impdats_smcfcs, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)




# Predicted probabilities CCA ---------------------------------------------





# Pool predicted probabilites ---------------------------------------------
# 
# tmat <- mstate::trans.comprisk(2, c("Rel", "NRM"))
# 
# dat_msprepped <- mstate::msprep(
#   time = c(NA, "ci_allo1", "ci_allo1"),
#   status = c(NA, "ev1", "ev2"), 
#   data = dat_mds,
#   trans = tmat,
#   keep = predictors
# ) 
# 
# dat_expanded <- mstate::expand.covs(
#   dat_msprepped, predictors, append = TRUE, longnames = F
# )
# 
# # Function generating reference patient...
# # First supply single-row df, but will all factor levels intact
# refpat <- na.omit(dat_mds)[1, ..predictors]
# refpat
# model.matrix(~ ., data = refpat)
# # Remove intercept, order columns by.., 
# # colnames = 
# 
# # Expand covariates - deal with variable names later with longnames
# dat_expanded <- mstate::expand.covs(
#   dat_msprepped, predictors, append = TRUE, longnames = F
# )
# 
# # Keep trans-specific vars (i.e. with .1 or .2, .K at the end)
# cond <- grep(x = colnames(dat_expanded), pattern = "*(\\.[1-9+]$)")
# 
# # Full mod
# rhs_mstate <- paste(
#   colnames(dat_expanded[cond]),
#   collapse = " + "
# )
# 
# form_mstate <- as.formula(paste0("Surv(time, status) ~ strata(trans) + ", rhs_mstate))
# 
# survival::coxph(form_mstate, dat_expanded)
# cca_rel


# Presenting results ------------------------------------------------------




