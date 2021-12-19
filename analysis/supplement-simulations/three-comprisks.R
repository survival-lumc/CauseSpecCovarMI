devtools::load_all()
baseline <- CauseSpecCovarMI::mds_shape_rates

# Generate X and Z
n <- 2000
covars <- gen_covars(n = n, X_type = "continuous", r = 0.5)

# Try the same different hazards (only different right?) with another one on top
ev1_pars <- list(
  "a1" = 1.5, 
  "h1_0" = 0.04,
  "b1" = .5, 
  "gamm1" = 1
)

ev2_pars <- list(
  "a2" = baseline[baseline$state == "NRM", "shape"], 
  "h2_0" = baseline[baseline$state == "NRM", "rate"], 
  "b2" = .5, 
  "gamm2" = .5
)

ev3_pars <- list(
  "a3" = 1, # constant
  "h3_0" = 0.05, 
  "b3" = 0.75, 
  "gamm3" = .5
)

# Generate the latent times
T1 <- rweibull_KM(
  n, 
  alph = ev1_pars$a1, 
  lam = ev1_pars$h1_0 * exp(drop(c(ev1_pars$b1, ev1_pars$gamm1) %*% t(as.matrix(covars))))
)
T2 <- rweibull_KM(
  n, 
  alph = ev2_pars$a2, 
  lam = ev2_pars$h2_0 * exp(drop(c(ev2_pars$b2, ev2_pars$gamm2) %*% t(as.matrix(covars))))
)
T3 <- rweibull_KM(
  n, 
  alph = ev3_pars$a3, 
  lam = ev3_pars$h3_0 * exp(drop(c(ev3_pars$b3, ev3_pars$gamm3) %*% t(as.matrix(covars))))
)

t_tilde <- pmin(T1, T2, T3)
eps <- data.table::fcase(
  t_tilde == T1, 1,
  t_tilde == T2, 2,
  t_tilde == T3, 3
)
admin_cens <- 10
time <- pmin(t_tilde, admin_cens)
delta <- as.numeric(time < admin_cens) * eps
table(delta)

dat <- cbind.data.frame(time, delta, covars)
dat <- dat[order(dat$time), ]

# Add cumulative hazards and interactions
dat$ev1 <- with(dat, as.numeric(delta == 1))
dat$ev2 <- with(dat, as.numeric(delta == 2))
dat$ev3 <- with(dat, as.numeric(delta == 3))

dat$H1 <- nelsaalen_timefixed(dat, "time", "ev1")
dat$H2 <- nelsaalen_timefixed(dat, "time", "ev2")
dat$H3 <- nelsaalen_timefixed(dat, "time", "ev3")

dat$H1_Z <- with(dat, H1 * Z)
dat$H2_Z <- with(dat, H2 * Z)
dat$H3_Z <- with(dat, H3 * Z)

head(dat)

# Check the good stuff
library(survival)
library(mstate)
coxph(Surv(time, delta == 1) ~ X + Z, data = dat) |>  coef()
coxph(Surv(time, delta == 2) ~ X + Z, data = dat) |>  coef()
coxph(Surv(time, delta == 3) ~ X + Z, data = dat) |>  coef()



# Prepare mstate code..
covs <- c("X", "Z")
tmat <- trans.comprisk(K = 3)
msdat <- msprep(
  time = c(NA, rep("time", 3)),
  status = with(dat, cbind(NA, delta == 1, delta == 2, delta == 3)),
  trans = tmat,
  data = dat,
  keep = covs
)
msdat_exp <- expand.covs(msdat, covs = covs, longnames = FALSE)
ms_form <- reformulate(
  response = "Surv(time, status)",
  termlabels = c(paste0("X.", seq_len(3)), paste0("Z.", seq_len(3)), "strata(trans)")
)
mod <- coxph(ms_form, data = msdat_exp)
summarise_ref_CCA(mod, analy = "full")




# We dont do prediction here right?


induce_missings(n = nrow(dat), dat = dat, mech = "MAR", p = 0.5, eta1 = -1)

# Just mar missing (dont vary hazards or coef effects, mech strength and prop missing yes?) 
# Or just show all coefficients facetted?

# Use ch12, ch12_int and smcfcs
meths <- mice::make.method(dat)
predmat <- mice::make.predictorMatrix(dat)
predmat[] <- 0L
