devtools::load_all()

library(survival)

#Read-in parameter
baseline <- readRDS(
  "./inst/testdata/MDS_shape_rates.rds"
)

shape_ev1 <- baseline[baseline$state == "REL", "shape"]
base_rate_ev1 <- baseline[baseline$state == "REL", "rate"]

# Below is 'different' condition 
#shape_ev1 <- 1.5
#base_rate_ev1 <- 0.04


# Parameter Weibull event 1
ev1_pars <- list(
  "a1" = shape_ev1, 
  "h1_0" = base_rate_ev1,
  "b1" = 1, # beta1 = 1 
  "gamm1" = 1
)

# Parameters Weibull event 2
ev2_pars <- list(
  "a2" = baseline[baseline$state == "NRM", "shape"], 
  "h2_0" = baseline[baseline$state == "NRM", "rate"], 
  "b2" = .5, 
  "gamm2" = .5
)

# Generate a dataset based on scenario
dat <- SimsCauseSpecCovarMiss::generate_dat(
  n = 25000,
  X_type = "continuous", 
  r = 0.5, 
  ev1_pars = ev1_pars,
  ev2_pars = ev2_pars, 
  rate_cens = baseline[baseline$state == "EFS", "rate"], 
  mech = "MAR", 
  p = 0.5,
  eta1 = 1
)

ev1_pars

coxph(Surv(t, eps == 1) ~ X + Z, data = dat)


dat$miss_ind
dat$X_indo <- with(dat, ifelse(is.na(X), 0, X))
mod_ind <- coxph(Surv(t, eps == 1) ~ X_indo + Z + miss_ind, data = dat)
mod_ind



# Generate a dataset based on scenario
dat <- SimsCauseSpecCovarMiss::generate_dat(
  n = 25000,
  X_type = "binary", 
  r = 0.5, 
  ev1_pars = ev1_pars,
  ev2_pars = ev2_pars, 
  rate_cens = baseline[baseline$state == "EFS", "rate"], 
  mech = "MCAR", 
  p = 0.5,
  eta1 = 1
)

ev1_pars

coxph(Surv(t, eps == 1) ~ X + Z, data = dat)

dat$X_copy <- ifelse(is.na(as.character(dat$X)), 2, as.character(dat$X))


dat$X_indo <- factor(dat$X_copy, levels = c(0, 1, 2))
mod_ind <- coxph(Surv(t, eps == 1) ~ X_indo + Z, data = dat)
mod_ind


# In practice; miss_ind for whether complete case or not?