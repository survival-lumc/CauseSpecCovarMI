##************************************************##
## Testing what may be responsible for small bias ##
##          in predicted probabilities            ##
##************************************************##


devtools::load_all()
library(tidyverse)
library(data.table)

set.seed(4328)

# Generate data -----------------------------------------------------------


# Read-in parameter
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
  n = 20000,
  X_type = "continuous", 
  r = 0.5, 
  ev1_pars = ev1_pars,
  ev2_pars = ev2_pars, 
  #rate_cens = baseline[baseline$state == "EFS", "rate"], 
  rate_cens = 0.01,
  mech = "MAR", 
  p = 0,
  eta1 = -1
)


# Checking the baseline hazards -------------------------------------------


tmat <- mstate::trans.comprisk(2, c("Rel", "NRM"))

# Fit model
mod <- setup_mstate(dat)

baseline_dat <- data.frame(X.1 = c(0, 0),
                           X.2 = c(0, 0), 
                           Z.1 = c(0, 0), 
                           Z.2 = c(0, 0), 
                           trans = c(1, 2), 
                           strata = c(1, 2))

msfit_baseline <- mstate::msfit(mod,
                                newdata = baseline_dat,
                                trans = tmat)

# We have no missing data in either eps or t, so happy days!
msfit_baseline$Haz %>% 
  data.table() %>% 
  .[, .N, by = trans]


test_probs <- function(b1,
                       b2,
                       Xval,
                       Zval,
                       ms_obj) {
  
  
  combs <- as.data.frame(expand.grid("b1" = b1, b2 = b2))
  
  hi <- purrr::map_dfr(1:nrow(combs), function(row) {
    
    test <- combs[row, ]
    b1 <- test$b1
    b2 <- test$b2
    
    dat <- get_state_probs(
      obj = ms_obj,
      combo = list("val_X" = Xval, 
                   "val_Z" = Zval),
      list("coefficients" = c(b1, 1, b2, .5))
    ) %>% 
      filter(time <= 5) %>% 
      slice(n()) %>% 
      mutate(times = 5) %>% 
      left_join(
        get_true_cuminc(ev1_pars, ev2_pars, 
                        list("val_X" = Xval,  
                             "val_Z" = Zval), 
                        times = 5), by = "times"
      ) %>% 
      select(times, 
             pstate1, true_pstate1,
             pstate2, true_pstate2,
             pstate3, true_pstate3) %>% 
      cbind.data.frame("b1" = b1,
                       "b2" = b2)
  })
  
  p <- data.table(hi) %>% 
    .[, ':=' (
      bias_EFS = pstate1 - true_pstate1,
      bias_REL = pstate2 - true_pstate2,
      bias_NRM = pstate3 - true_pstate3,
      b1 = factor(b1),
      b2 = factor(b2),
      inter = interaction(b1, b2)
    )] %>% 
    melt.data.table(
      id.vars = c("inter", "b1", "b2"), 
      measure.vars = c("bias_EFS", "bias_REL", "bias_NRM"), 
      variable.name = "state",
      value.name = "bias", 
    ) %>% 
    ggplot(aes(b1, b2, fill = abs(bias))) +
    geom_tile() +
    facet_wrap(. ~ state) +
    scale_fill_gradient(low="white", high="blue") + 
    coord_cartesian(expand = 0) +
    ggtitle("Data-generating b1 = 1, and b2 = 0.5; pred at 5y")
  
  return(p)
}

# Check baseline hazards
test_probs(b1 = seq(0.8, 1.2, by = 0.1),
           b2 = seq(0.3, 0.7, by = 0.1),
           Xval = 0,
           Zval = 0,
           ms_obj = msfit_baseline)

# Check a combo
test_probs(b1 = seq(0.8, 1.2, by = 0.1),
           b2 = seq(0.3, 0.7, by = 0.1),
           Xval = .5,
           Zval = .5,
           ms_obj = msfit_baseline)




library(CFC)
library(survival)

# Models without CFC

#
# https://stackoverflow.com/questions/32856440/simulate-competing-risk-data

# This is not supported
survreg(Surv(t, eps) ~ 1, data = dat, 
        dist = "weibull")


# If you remove artificial censoring this is all correct
mod_cens <- survreg(Surv(t, eps == 0) ~ 1, data = dat, dist = "exponential")
mod_rel <- survreg(Surv(t, eps == 1) ~ X + Z, data = dat, dist = "weibull")
mod_nrm <- survreg(Surv(t, eps == 2) ~ X + Z, data = dat, dist = "weibull")

shape_rate <- function(mod) {
  
  mod_coefs <- mod$coefficients
  shape <- 1 / mod$scale
  rate <- exp(mod_coefs[1])^(-shape)

  cat(shape, rate)
}

shape_rate(mod_cens)
shape_rate(mod_rel)
c(ev1_pars[[1]], ev1_pars[[2]])
shape_rate(mod_nrm)
c(ev2_pars[[1]], ev2_pars[[2]])



# Try with cfc
formul <- Surv(t, eps) ~ X
hi <- cfc.survreg(formula = formul, data = dat, newdata = dat, dist = "weibull") 

?cfc.survreg

dat$eps <- as.integer(dat$eps) - 1

data(bmt)
formul <- Surv(time, cause) ~ platelet + age + tcell
ret <- cfc.survreg(formul, bmt[1:300, ], bmt[-(1:300), ]
                   , Nmax = 300, rel.tol = 1e-3)

bmt$cause

# Probability missing demo ------------------------------------------------


# The following holds for 50% missing - because Z is standard normal
# intercept (eta0) will always be approx 0


etas <- c(0, -.25, -0.5, -.75, -1, -2)

Z <- rnorm(1000)

prob_plots <- purrr::map_dfr(
  etas,
  ~ cbind.data.frame(
    "Z" = Z, 
    "prob" = plogis(.x * Z),
    "eta1" = .x
  )
)

prob_plots %>% 
  dplyr::mutate(eta1 = factor(eta1)) %>% 
  ggplot(aes(Z, prob, col = eta1, group = eta1)) +
  geom_line(size = 1) +
  theme_bw() +
  coord_cartesian(expand = 0.1)


# Facetted plot

test <- purrr::map_dfr(
  etas,
  ~ generate_dat(
    n = 500,
    X_type = "continuous", 
    r = 0.5, 
    ev1_pars = ev1_pars,
    ev2_pars = ev2_pars, 
    rate_cens = baseline[baseline$state == "EFS", "rate"], 
    mech = "MAR_GEN", 
    p = 0.5,
    eta1 = .x
  ) %>% 
    dplyr::mutate("eta1" = .x)
)

test %>% 
  dplyr::mutate(
    miss_ind = factor(miss_ind),
    eta1 = factor(eta1, levels = as.character(etas),
                  labels = paste0("eta1 = ", as.character(etas)))
  ) %>% 
  ggplot(aes(t, X_orig, col = miss_ind)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_colour_manual("Missing indicator",
                      labels = c("Observed", "Missing"),
                      values = c(9, 6)) +
  xlab("t") +
  facet_wrap(. ~ eta1) +
  theme_bw()



# Example with rnorm around 5
Z <- rnorm(1000, mean = 5)
plot(Z, plogis(-5 + Z))

mean(plogis(Z))
