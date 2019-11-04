# Load libraries
library(smcfcs)
library(mitools)

source("data_generation/dat_generation_W.R")



# Parameters 
true_betas <- c(0.5, -0.5) # beta (X1) and gamma (X2) in data generating model
mech <- "MAR" # missingness mechanism
m <- c(1, 2, 5) # number of imputations to compare
method <- "norm" # imputation method, here bayesian linear regression
prob <- 0.3 # percentage missingness in X!

# Generate 
dat <- dat_gener_W(N = 200, 
                   X_type = "contin",
                   mus = c(0, 0), 
                   covmat = matrix(c(1, 0.25, 
                                     0.25, 1), nrow = 2), # make into correlation mat
                   mech = mech, 
                   pars_MAR = c(1, 1),
                   p = prob, 
                   cause2 = "weib", 
                   vals_t1 = c(0.3, 1, true_betas), 
                   vals_t2 = c(1.7, 0.5, -0.5, 0.5))

# Attempt
imps <- quiet(smcfcs(originaldata = dat, 
                     smtype = "compet", 
                     smformula = c("Surv(t, eps == 1) ~ X1 + X2",
                                   "Surv(t, eps == 2) ~ X1 + X2"), 
                     method = c("norm", rep("", 10)), 
                     m = m[length(m)]))



