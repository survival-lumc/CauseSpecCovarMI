##***************************##
## Generating synthetic data ##
##***************************##


# (For possible sharing..) - based on synthpop
options(contrasts = rep("contr.treatment", 2)) 

# Read-in (without admin censoring)
dat_mds <- fst::read_fst("data/dat-mds.fst") %>% 
  setDT()


# Prepare synthpop --------------------------------------------------------


dat_mds[, ':=' (
  srv_s_allo1 = NULL,
  srv_allo1 = NULL
)]

x <- dat_mds[, !c("ci_allo1", "ci_s_allo1")]
test_syn <- synthpop::syn(x)
xp <- test_syn$syn


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


synth_covar <- synthpop::syn(
  data = dat_mds,
  method = meths, 
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
  imp_type = "smcfcs",
  cont_method = "norm"
)

# Make formula
smform_smcfcs <- c(
  Reduce(paste, deparse(form_rel)),
  Reduce(paste, deparse(form_nrm))
)

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
  )] %>% 
  
  # Gen new event times
  .[, ':=' (
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


# Save synthetic data -----------------------------------------------------


fst::write_fst(x = dat_mds_synth$syn, path = "analysis/data/dat-mds-synth.fst")
