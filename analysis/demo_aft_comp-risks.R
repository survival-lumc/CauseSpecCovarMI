

library(survival)
devtools::load_all()

set.seed(1984)

# Sample size
n <- 25000

# Covariates
covars <- MASS::mvrnorm(
  n = n, 
  mu = c(1, 3), 
  Sigma = matrix(
    c(1, 0.5,
      0.5, 1), 
    nrow = 2,
    ncol = 2
  )
) %>% 
  magrittr::set_colnames(c("X", "Z")) %>% 
  as.data.frame()


# Events 1
shape_ev1 <- 0.63
base_ev1 <- 0.25
rate_ev1 <- with(covars, base_ev1 * exp(0.5 * X + -0.5 * Z))
t1 <- rweibull_KM(n = n, alph = shape_ev1, lam = rate_ev1)

# Event 2
shape_ev2 <- 1.44
base_ev2 <- 0.47
rate_ev2 <- with(covars, base_ev2 * exp(0.25 * X + -1 * Z))
t2 <- rweibull_KM(n = n, alph = shape_ev2, lam = rate_ev2)

# Censoring
rate_cens <- 0.13
cens <- rexp(n = n, rate = rate_cens)

# Generate dat
delta <- ifelse(t1 < t2, 1, 2)
t_tilde <- pmin(t1, t2)
t <- pmin(t_tilde, cens)
eps <- ifelse(cens < t_tilde, 0, delta)

# bind
dat <- cbind.data.frame(
  "t" = t,
  "eps" = eps,
  "X" = covars$X,
  "Z" = covars$Z
)

table(dat$eps)


# If you remove artificial censoring in generate_dat() this is all correct
mod_cens <- survreg(Surv(t, eps == 0) ~ 1, data = dat, dist = "exponential")
mod_rel <- survreg(Surv(t, eps == 1) ~ X + Z, data = dat, dist = "weibull")
mod_nrm <- survreg(Surv(t, eps == 2) ~ X + Z, data = dat, dist = "weibull")

# Without artif. censoring, this holds
shape_rate <- function(mod) {
  
  mod_coefs <- mod$coefficients
  shape <- 1 / mod$scale
  rate <- exp(mod_coefs[1])^(-shape)
  
  c(shape, rate)
}

# Check censoring
shape_rate(mod_cens)
rate_cens

# Ev1
shape_rate(mod_rel)
c(shape_ev1, base_ev1)

# Ev2
shape_rate(mod_nrm)
c(shape_ev2, base_ev2)

# Get prop haz coeffs
mod_rel$coefficients * (-shape_rate(mod_rel)[1])
mod_nrm$coefficients * (-shape_rate(mod_nrm)[1])
