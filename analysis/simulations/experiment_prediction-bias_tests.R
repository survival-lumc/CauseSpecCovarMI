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
    
devtools::load_all()
set.seed(4328)

# Generate a dataset based on scenario
dat <- SimsCauseSpecCovarMiss::generate_dat(
  n = 5000,
  X_type = "binary", 
  r = 0.5, 
  ev1_pars = ev1_pars,
  ev2_pars = ev2_pars, 
  rate_cens = baseline[baseline$state == "EFS", "rate"], 
  #rate_cens = 0.026,
  mech = "MAR_GEN", 
  p = 0.5,
  eta1 = -1
)

setup_mstate(dat) 
dat$X

coxph(Surv(t, eps == 1) ~ X_orig + Z, data = dat)

# -0.4810937
cor(dat$Z, dat$t)

summary(
  survival::coxph(Surv(t, eps == 1) ~ X + Z, data = dat)
)$coefficients[1,3]



lol <- replicate(n = 250, expr = {
  dat <- SimsCauseSpecCovarMiss::generate_dat(
    n = 500,
    X_type = "continuous", 
    r = 0.5, 
    ev1_pars = ev1_pars,
    ev2_pars = ev2_pars, 
    rate_cens = baseline[baseline$state == "EFS", "rate"], 
    #rate_cens = 0.026,
    mech = "MAR_GEN", 
    p = 0.5,
    eta1 = -1
  )
  
  
  
  summary(
    survival::coxph(Surv(t, eps == 1) ~ X + Z, data = dat)
  )$coefficients[1,3]
  
})

mean(lol)

hist(lol)

# Checking the baseline hazards -------------------------------------------


tmat <- mstate::trans.comprisk(2, c("Rel", "NRM"))

tmat <- mstate::trans.comprisk(2, c("Rel", "NRM"))
covs <- c("X", "Z")
covs <- c("X_fac", "Z")


dat$X_fac <- with(
  dat, cut(X_orig, breaks = c(-Inf, -1, 1.5, Inf), 
           labels = c("low", "med", "high"))
)

# Long format
dat_msprepped <- mstate::msprep(time = c(NA, "t", "t"),
                                status = c(NA, "ev1", "ev2"), 
                                data = dat,
                                trans = tmat,
                                keep = covs) 

# Expand covariates
dat_expanded <- mstate::expand.covs(dat_msprepped, covs,
                                    append = TRUE, longnames = T)


cox_long <- survival::coxph(Surv(time, status) ~ 
                              X.1 + Z.1 + # Trans == 1
                              X.2 + Z.2 + # Trans == 2
                              strata(trans), # Separate baseline hazards
                            data = dat_expanded)


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

library(survival)


# If you remove artificial censoring in generate_dat() this is all correct
mod_cens <- survreg(Surv(t, eps == 0) ~ 1, data = dat, dist = "exponential")
mod_rel <- survreg(Surv(t, eps == 1) ~ X + Z, data = dat, dist = "weibull")
mod_nrm <- survreg(Surv(t, eps == 2) ~ X + Z, data = dat, dist = "weibull")

# Without artif. censoring, this holds
shape_rate <- function(mod) {
  
  mod_coefs <- mod$coefficients
  shape <- 1 / mod$scale
  rate <- exp(mod_coefs[1])^(-shape)

  cat(shape, rate)
}

shape_rate(mod_cens)#0.026
shape_rate(mod_rel)
c(ev1_pars[[1]], ev1_pars[[2]])
shape_rate(mod_nrm)
c(ev2_pars[[1]], ev2_pars[[2]])

library(prodlim)
library(riskRegression)

fit1 <- CSC(list(Hist(t, eps)) ~ X + Z, data=dat)

library("flexsurv")

tmat <- mstate::trans.comprisk(2, c("Rel", "NRM"))
covs <- c("X", "Z")

# Long format
dat_msprepped <- mstate::msprep(time = c(NA, "t", "t"),
                                status = c(NA, "ev1", "ev2"), 
                                data = dat,
                                trans = tmat,
                                keep = covs) 

# Expand covariates
dat_expanded <- mstate::expand.covs(dat_msprepped, covs,
                                    append = TRUE, longnames = F)


test_flex <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), 
                         data = dat_msprepped,
                         dist = "weibullPH")

flexsurvreg(Surv(time, status) ~ X + Z, 
            subset = (trans == 2),
            data = dat_msprepped,
            dist = "weibullPH")



test_flex

coxph(Surv(t, eps) ~ X + Z, data = dat, id = 1:nrow(dat))

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

# Try with onko kiss dara
#install.packages("../../software-other/compeir_1.0.tar.gz", repos = NULL, type="source")
library("compeir")
data("okiss")

# Mod 

#
library(tidyverse)
library(survival)

dato <- okiss %>% 
  mutate(eps = case_when(
    status == 11 ~ 0,
    status == 1 ~ 1,
    status == 2 ~ 2,
    status == 7 ~ 2,
  ))

# mod BSI
survreg(Surv(time, eps == 0) ~ 1,
        data = dato, dist = "lognormal")

survreg(Surv(time, eps == 1) ~ 1,
        data = dato, dist = "lognormal")

survreg(Surv(time, eps == 2) ~ 1,
        data = dato, dist = "lognormal")


survreg(Surv(time, eps == 1) ~ sex + allo,
        data = dato, dist = "lognormal")

survreg(Surv(time, eps == 2) ~ sex + allo,
        data = dato, dist = "lognormal")


# Probability missing demo ------------------------------------------------


# The following holds for 50% missing - because Z is standard normal
# intercept (eta0) will always be approx 0


etas <- c(0, -.25, -0.5, -.75, -1, -2)

Z <- seq(-3, 3, by = 0.01)

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
  geom_line(size = 1.5) +
  coord_cartesian(expand = 0, ylim = c(0, 1)) +
  ylab(bquote(P(R[X] ~ "= 1"*"|"*Z))) +
  guides(
    col = guide_legend(title = expression(eta[1]))
  ) +
  scale_color_brewer(palette = "Dark2") 


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
    mech = "MAR", 
    p = 0.5,
    eta1 = .x
  ) %>% 
    dplyr::mutate("eta1" = .x)
)

expression(paste(eta, as.character(etas)))


test %>% 
  dplyr::mutate(
    miss_ind = factor(miss_ind),
    eta1 = factor(eta1, levels = as.character(etas),
                  labels = paste0("eta1 = ", as.character(etas)))
  ) %>% 
  ggplot(aes(Z, X_orig, col = miss_ind)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_colour_manual("Missingness indicator",
                      labels = c("Observed", "Missing"),
                      values = c(9, 6)) +
  xlab("Z") +
  ylab("X") +
  facet_wrap(. ~ eta1) +
  theme(legend.position = "top")



test %>% 
  dplyr::mutate(
    miss_ind = factor(miss_ind),
    eta1 = factor(eta1, levels = as.character(etas),
                  labels = paste0("eta1 = ", as.character(etas)))
  ) %>% 
  filter(eta1 %in% c("eta1 = 0", "eta1 = -0.5", "eta1 = -2")) %>% 
  ggplot(aes(Z, X_orig, col = miss_ind)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_colour_manual("",
                      labels = c("Observed", "Missing"),
                      values = c(9, 6)) +
  xlab("Z") +
  ylab("X") +
  facet_wrap(. ~ eta1, nrow = 3) +
  theme_bw(base_size = 30)+
  theme(legend.position = "top") 


ggsave("MAR_X-Z.svg", dpi = "retina", units = "in",
       width = 8, height = 10)


# Example with rnorm around 5
Z <- rnorm(1000, mean = 5)
plot(Z, plogis(-5 + Z))

mean(plogis(Z))


# Lets do same procedure with t
etas <- c(0, -.25, -0.5, -.75, -1, -2)

t <- seq(0.01, max(dat$t), by = 0.1)

prob_plots <- purrr::map_dfr(
  etas,
  ~ cbind.data.frame(
    "t" = dat$t, 
    "prob" = plogis(.x * scale(log(dat$t))),
    "eta1" = .x
  )
)

prob_plots %>% 
  dplyr::mutate(eta1 = factor(eta1)) %>% 
  ggplot(aes(t, prob, col = eta1, group = eta1)) +
  geom_line(size = 1.5) +
  coord_cartesian(expand = 0, ylim = c(0,1)) +
  ylab(bquote(P(R[X] ~ "= 1"*"|"*T[stand]))) +
  guides(
    col = guide_legend(title = expression(eta[1]))
  ) +
  xlab("T") +
  scale_color_brewer(palette = "Dark2") 


prob_plots %>% 
  dplyr::mutate(eta1 = factor(eta1)) %>% 
  ggplot(aes(scale(log(t)), prob, col = eta1, group = eta1)) +
  geom_line(size = 1.5) +
  theme_bw() +
  coord_cartesian(expand = 0.1) +
  ylab(bquote(P(R[X] ~ "= 1"*"|"*T[stand]))) +
  guides(
    col = guide_legend(title = expression(eta[1]))
  ) +
  xlab(bquote(T[stand])) + 
  scale_color_brewer(palette = "Dark2")


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
  xlab("T") +
  ylab("X") +
  facet_wrap(. ~ eta1) +
  theme(legend.position = "top")



# Plot illust haz ---------------------------------------------------------



# Plots of cumulative incidences and hazards
theme_set(theme_bw(base_size = 24))
library(tidyverse)


combo = data.frame(
  "val_X" = 0, 
  "val_Z" = 0
)                  

ev1_pars = list(
  "a1" = 0.58, 
  "h1_0" = 0.19,
  "b1" = 0, 
  "gamm1" = 0
)    

ev2_pars = list(
  "a2" = 0.53, 
  "h2_0" = 0.21, 
  "b2" = 0, 
  "gamm2" = 0
)


cuminc_similar <- get_true_cuminc(ev1_pars, ev2_pars, 
                combo, 
                times = seq(0, 10, by = 0.25)) %>% 
  pivot_longer(true_pstate2:true_pstate1,
               names_to = "state",
               values_to = "prob") %>% 
  ggplot(aes(times, prob, col = state,
             linetype = state)) +
  geom_line(size = 1.5, alpha = .75) +
  scale_colour_manual(
    "State", 
    values = c(1, 2, 3),
    labels = c("EFS", "REL", "NRM")
  ) +
  scale_linetype_manual("State", values = c("dotted", "solid","solid"),
                        labels = c("EFS", "REL", "NRM")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  
  # Title and labs
  xlab("Time (years)") + 
  ylab("Probability")  +
  ggtitle("haz_shape = 'similar'")


t <- seq(0.01, 10, by = 0.1)

lam1 <- ev1_pars$h1_0
lam2 <- ev2_pars$h2_0

haz1 <- haz_weib(alph = ev1_pars$a1, lam = lam1, t = t) 
haz2 <- haz_weib(alph = ev2_pars$a2, lam = lam2, t = t) 
haz_cens <- haz_weib(alph = 1, lam = 0.14, t = t)


haz_similar <- cbind.data.frame(
  "t" = t,
  "haz1" = haz1,
  "haz2" = haz2,
  "haz3" = haz_cens
) %>% 
  gather("event", "hazard", .data$haz1, .data$haz2, .data$haz3) %>% 
  ggplot(aes(.data$t, .data$hazard, col = event,
             linetype = event)) +
  geom_line(size = 1.5, alpha = 0.7) +
  xlab("Time") +
  ylab("Hazard") +
  scale_color_manual(
    "State", 
    values = c(2, 3, 1),
    labels = c("REL", "NRM", "EFS")
  ) +
  scale_linetype_manual("State", values = c("solid", "solid","dotted"),
                        labels = c("REL", "NRM", "EFS")) + 
  theme(legend.position = "bottom")





ev1_pars = list(
  "a1" = 1.5, 
  "h1_0" = 0.04,
  "b1" = 0, 
  "gamm1" = 0
)    

ev2_pars = list(
  "a2" = 0.53, 
  "h2_0" = 0.21, 
  "b2" = 0, 
  "gamm2" = 0
)


cuminc_diff <-get_true_cuminc(ev1_pars, ev2_pars, 
                combo, 
                times = seq(0, 10, by = 0.25)) %>% 
  pivot_longer(true_pstate2:true_pstate1,
               names_to = "state",
               values_to = "prob") %>% 
  ggplot(aes(times, prob, col = state,
             linetype = state)) +
  geom_line(size = 1.5, alpha = .75) +
  scale_colour_manual(
    "State", 
    values = c(1, 2, 3),
    labels = c("EFS", "REL", "NRM")
  ) +
  scale_linetype_manual("State", values = c("dotted", "solid","solid"),
                        labels = c("EFS", "REL", "NRM")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  
  # Title and labs
  xlab("Time (years)") + 
  ylab("Probability") +
  ggtitle("haz_shape = 'different'")

t <- seq(0.01, 10, by = 0.1)

lam1 <- ev1_pars$h1_0
lam2 <- ev2_pars$h2_0

haz1 <- haz_weib(alph = ev1_pars$a1, lam = lam1, t = t) 
haz2 <- haz_weib(alph = ev2_pars$a2, lam = lam2, t = t) 
haz_cens <- haz_weib(alph = 1, lam = 0.14, t = t)

haz_diff <- cbind.data.frame(
  "t" = t,
  "haz1" = haz1,
  "haz2" = haz2,
  "haz3" = haz_cens
) %>% 
  gather("event", "hazard", .data$haz1, .data$haz2, .data$haz3) %>% 
  ggplot(aes(.data$t, .data$hazard, col = event,
             linetype = event)) +
  geom_line(size = 1.5, alpha = 0.7) +
  xlab("Time") +
  ylab("Hazard") +
  scale_color_manual(
    "State", 
    values = c(2, 3, 1),
    labels = c("REL", "NRM", "EFS")
  ) +
  scale_linetype_manual("State", values = c("solid", "solid","dotted"),
                        labels = c("REL", "NRM", "EFS")) + 
  theme(legend.position = "bottom")



# library ggpubr
library(ggpubr)

ggarrange(cuminc_similar, cuminc_diff,
          haz_similar, haz_diff, ncol = 2, nrow = 2, 
          common.legend = T,
          legend = "bottom")


p <- ggarrange(cuminc_similar,
          cuminc_diff, nrow = 2,
          legend = "bottom",
          common.legend = T) 

p

ggsave("cumincs.svg", dpi = "retina", units = "in",
       width = 8, height = 9)
