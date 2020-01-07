
#
library(JointAI)
library(MASS)
library(tidyverse)
library(survival)
library(mice)
library(mstate)
library(cmprsk)

source("shark/data_generation/dat_generation_MRR.R")

dat <- dat_gener_MRR(N = 500, 
                   X_type = "contin",
                   mus = c(0, 0), 
                   covmat = matrix(c(1, 0.25, 
                                     0.25, 1), nrow = 2), # make into correlation mat
                   mech = "MAR", 
                   #pars_MAR = c(1, 1),
                   p = .2, 
                   cause2 = "weib", 
                   vals_t1 = c(0.3, 1, c(0.5, -0.5)), 
                   vals_t2 = c(1.7, 0.5, -0.5, 0.5)) 

dat <- dat_gener_W(N = 500, 
                     X_type = "contin",
                     mus = c(0, 0), 
                     covmat = matrix(c(1, 0.25, 
                                       0.25, 1), nrow = 2), # make into correlation mat
                     mech = "MAR", 
                     pars_MAR = c(1, 1),
                     p = .2, 
                     cause2 = "weib", 
                     vals_t1 = c(0.3, 1, c(0.5, -0.5)), 
                     vals_t2 = c(1.7, 0.5, -0.5, 0.5)) 

dat$ind <- factor(as.numeric(dat$eps == 1), levels = c(0, 1))

indo_X2 <- sample(T:F, size = 500, T, prob = c(0.2, 0.8))
dat$X2 <- ifelse(indo_X2, NA, dat$X2)

# With joint ai
mod <- coxph_imp(Surv(t, ind) ~ X1 + X2, data = dat, n.chains = 3,
                 n.iter = 100, keep_model = T,
                 monitor_params = c("analysis_main" = TRUE,
                                    "imps" = TRUE,
                                    "other" = "log.h0.T")) #log.h0.t

mod$model

mod$model$data()

list_models(mod)

traceplot(mod)

# JAGS model
mod$MCMC

mod$monitor_params

colnames(mod$MCMC[[1]])

mod$MCMC[[1]]

# https://cran.r-project.org/web/packages/JointAI/vignettes/TheoreticalBackground.html

#as.data.frame(mod$MCMC[[1]]) %>% 
bind_rows(as.data.frame(mod$MCMC[[1]]),
          as.data.frame(mod$MCMC[[2]]),
          as.data.frame(mod$MCMC[[3]])) %>% 
  select("log.h0.T[1]":"log.h0.T[500]") %>% 
  gather(id, h0) %>%
  #mutate(samp = rep(1:100, 500)) %>%  # index for posterior draw
  mutate(samp = rep(1:300, 500)) %>%  # index for posterior draw    
  #https://stackoverflow.com/questions/2403122/regular-expression-to-extract-text-between-square-brackets
  # Answer i slets breaks it down
  mutate(id = as.numeric(str_extract(id, "(?<=\\[).*?(?=\\])"))) %>% 
  group_by(samp) %>% 
  mutate(t = dat$t) %>% 
  ungroup() %>% 
  
  # Use exp(h0)
  ggplot(aes(t, h0)) + geom_line(aes(col = factor(samp))) +
  theme(legend.position = "none")


# 
bind_rows(as.data.frame(mod$MCMC[[1]]),
         as.data.frame(mod$MCMC[[2]]),
         as.data.frame(mod$MCMC[[3]])) %>% 
  mutate(iter = 1:n()) %>% 
  ggplot(aes(iter, `Xc[125,3]`)) + geom_line()



# Get imped datasets
test <- get_MIdat(mod, m = 4, include = F, start = 1, )
test



# Find imputations
as.data.frame(mod$MCMC[[1]])



ncol(mod$MCMC[[1]])

# use minspace
test$Imputation_

# Full data
cox1 <- coxph(Surv(t, eps == 1) ~ X1_orig + X2, data = dat)

dat$event <- factor(dat$eps, levels = c(0, 1, 2),
                    labels = c("censor", "REL", "NRM"))
dat$id <- 1:nrow(dat)

fito <- coxph(Surv(t, event) ~ X1_orig + X2, data=dat, id = id)

# 
dummy <- data.frame("X1_orig" = 0, "X2" = 0)
predicto <- survfit(fito, newdata = dummy)
plot(predicto)

predicto_cox1 <- survfit(cox1, newdata = dummy)
plot(predicto_cox1, conf.int = F, fun = function(x) 1 - x)

# Using mstate
tmat <- trans.comprisk(2, names = c("cens", "REL", "NRM"))
tmat

dat$stat1 <- as.numeric(dat$eps == 1)
dat$stat2 <- as.numeric(dat$eps == 2)

dat_long <- msprep(time = c(NA, "t", "t"),
                   status = c(NA, "stat1", "stat2"), 
                   data = dat, 
                   keep = c("X1_orig", "X2"), trans = tmat)

events(dat_long)

dat_long <- expand.covs(dat_long, c("X1_orig", "X2"))

coxph(Surv(time, status) ~ X1_orig + X2, dat = dat_long,
      subset = (trans == 2))

c1 <- coxph(Surv(time, status) ~ X1_orig.1 + X2.1 + 
              X1_orig.2 + X2.2 + strata(trans), 
            data = dat_long,
            method = "breslow")

c1

WW <- data.frame(X1_orig.1 = c(0, 0),
                 X2.1 = c(0, 0),
                 X1_orig.2 = c(0, 0),
                 X2.2 = c(0, 0), 
                 trans = c(1, 2), strata = c(1, 2))

msf.WW <- msfit(c1, WW, trans = tmat)

pt.WW <- probtrans(msf.WW, 0)[[1]]


