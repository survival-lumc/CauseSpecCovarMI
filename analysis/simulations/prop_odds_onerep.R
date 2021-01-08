
# Set a seed
set.seed(1984)
options(contrasts = rep("contr.treatment", 2)) 
devtools::load_all()


# Set params
n <- 20000
Z <- rnorm(n = n, mean = 0, sd = 1)

# Generate X
X_latent <- 1 * Z + rlogis(n = n, location = 0, scale = 1)
plot(X_latent, Z)
cor(X_latent, Z)

# Define break points and get X - based of hctci risks (disprop larger lower half)
break_points <- c(-Inf, 0, 1.5, Inf)
X <- cut(X_latent, breaks = break_points, ordered_result = T, 
         labels = c("low", "intermediate", "high"))


round(prop.table(table(X)), 2)
round(prop.table(table(dat_mds_reg$cytog_threecat)), 2)

dat_covars <- cbind.data.frame(X, Z)


gen_cmprsk_times(
  n = n,
  dat = dat_covars,
  ev1_pars = ev1_pars,
  ev2_pars = ev2_pars, 
  rate_cens = 0.0001 
)


# Generate competing risks data
baseline <- readRDS(
  "./inst/testdata/MDS_shape_rates.rds"
)

shape_ev1 <- 1.5
base_rate_ev1 <- 0.04
scenario <- list("beta1" = 1)

# Parameter Weibull event 1
ev1_pars <- list(
  "a1" = shape_ev1, 
  "h1_0" = base_rate_ev1,
  "b1" = c(0.5, scenario$beta1), 
  "gamm1" = 1
)

# Parameters Weibull event 2
ev2_pars <- list(
  "a2" = baseline[baseline$state == "NRM", "shape"], 
  "h2_0" = baseline[baseline$state == "NRM", "rate"], 
  "b2" = c(.5, .5), 
  "gamm2" = .5
)


#
dat <- SimsCauseSpecCovarMiss::generate_dat(
  n = 2000,
  X_type = "ordcat", 
  ev1_pars = ev1_pars,
  ev2_pars = ev2_pars, 
  rate_cens = baseline[baseline$state == "EFS", "rate"], 
  mech = "MAR", 
  p = 0.5,
  eta1 = -2, 
  mod_type = "latent"
)




setup_mstate(dat %>%  mutate(X = X_orig))
setup_mstate(dat)

unlist(ev1_pars)
coxph(Surv(t, eps == 1) ~ X_orig + Z, data = dat)

unlist(ev2_pars)
coxph(Surv(t, eps == 2) ~ X_orig + Z, data = dat)

dat %>% 
  ggplot(aes(X_orig, Z, col = factor(miss_ind))) +
  geom_jitter(width = 0.2, alpha = 0.75)

dat %>% 
  ggplot(aes(X_orig, t, col = factor(miss_ind))) +
  geom_jitter(width = 0.2, alpha = 0.75)


# 
mod_mat <- model.matrix(~ X + Z, contrasts.arg = list(X = "contr.treatment"))
lam1 <- ev1_pars$h1_0 * exp(mod_mat %*% c(0, ev1_pars$b1, ev1_pars$gamm1))
lam2 <- ev2_pars$h2_0 * exp(mod_mat %*% c(0, ev2_pars$b2, ev2_pars$gamm2))

t1 <- rweibull_KM(n = n, alph = ev1_pars$a1, lam = lam1)
t2 <- rweibull_KM(n = n, alph = ev2_pars$a2, lam = lam2)

# Take minimum of two, and compute indicator
t <- pmin(t1, t2)
eps <- ifelse(t1 < t2, 1, 2)

cens <- stats::rexp(n = n, rate = baseline[baseline$state == "EFS", "rate"])
eps <- ifelse(cens < t, 0, eps)
t <- pmin(cens, t)

# Also artificially censor at 10 years
eps <- ifelse(t >= 10, 0, eps)
t <- pmin(t, 10)

# Bind - keep latent to use for MNAR
dat <- cbind.data.frame(X, X_latent, Z, eps, t)

blop <- induce_missings(n, dat, p = 0.5, mech = "MNAR", eta1 = -1)

prop.table(table(blop$X))
prop.table(table(blop$X_miss))


blop %>% 
  ggplot(aes(X, Z, col = factor(miss_ind))) +
  geom_jitter(width = 0.2)

ev1_pars

coxph(Surv(t, eps == 2) ~ factor(X_miss) + Z, data = blop)
