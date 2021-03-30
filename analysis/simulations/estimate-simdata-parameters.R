##*****************************************##
## Estimating AFT parameters from MDS data ##
##*****************************************##


# Read-in
dat_mds <- fst::read_fst(path = "data/dat-mds_admin-cens.fst") %>% 
  data.table::setDT()


# Relapse -----------------------------------------------------------------


aft_rel <- survival::survreg(
  formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ 1, 
  data = dat_mds, 
  dist = "weibull"
)

shape_rel <- 1 / aft_rel$scale
rate_rel <- exp(aft_rel$coefficients)^(-shape_rel)


# NRM ---------------------------------------------------------------------


aft_nrm <- survival::survreg(
  formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ 1, 
  data = dat_mds, 
  dist = "weibull"
)

coefs_nrm <- aft_nrm$coefficients
shape_nrm <- 1 / aft_nrm$scale
rate_nrm <- exp(aft_nrm$coefficients)^(-shape_nrm)


# Censoring ---------------------------------------------------------------


# Censoring is event, and those with admin censoring are set to event
# just after 10 years
dat_mds[, srv_s_rev := ifelse(srv_allo1 == 10, 1, srv_s_allo1)]
dat_mds[, srv_allo1_rev := ifelse(srv_allo1 == 10, srv_allo1 + 0.01, srv_allo1)]


# Exponential AFT this time
aft_cens <- survival::survreg(
  formula = Surv(srv_allo1_rev, srv_s_rev == 0) ~ 1, 
  data = dat_mds, 
  dist = "exponential"
)

rate_cens <- exp(aft_cens$coefficients)^(-1)


# Save parameters in df ---------------------------------------------------

mds_shape_rates <- cbind.data.frame(
  "state" = c("REL", "NRM", "EFS"),
  "shape" = c(shape_rel, shape_nrm, 1),
  "rate" = c(rate_rel, rate_nrm, rate_cens)
)

mds_shape_rates

# Write into fst (write original rds one here)
#fst::write_fst(x = mds_shape_rates, path = "analysis/data/mds-shape-rates.fst")
