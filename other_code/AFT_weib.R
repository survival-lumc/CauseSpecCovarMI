# AFT demo (have slides at the ready)

source("data_generation/dat_generation_MRR_KMweib.R")

# Generate single time to event
t <- my_rweibull(10000, alph = 1.25, lam = 3)
mod <- survreg(Surv(t) ~ 1, dist = "weibull") # fit a baseline model

coefs <- mod$coefficients # same as mod$icoef

# See help file of survreg
shap <- 1 / mod$scale 
shap # this is the shape 
exp(coefs[1])^(-shap) # This is the rate, in KM param need to take power of -shap


# Try with simulated data
dat <- dat_gener_MRR_KM(N = 100000,
                        X_type = "contin",
                        mus = c(0, 0), 
                        covmat = matrix(c(1, 0.25, 
                                          0.25, 1), nrow = 2), 
                        mech = "MCAR", 
                        p = 0, 
                        cause2 = "weib", 
                        vals_t1 = c(0.3, 5, 1, -0.5), 
                        vals_t2 = c(1.7, 0.5, -0.5, 0.5))

# Cause specific models
coxph(Surv(t, eps == 1) ~ X1 + X2, data = dat)

# Survreg with both parameters
mod1 <- survreg(Surv(t, eps == 1) ~  X1 + X2, data = dat, dist = "weibull")

# Two objects:
mod1$icoef # parameters of marginal modelsurvreg(Surv(t, eps == 1) ~  1, data = dat, dist = "weibull")
mod1$coefficients # model coefficients


# Lets get back coefficients
shap <- 1 / mod1$scale
shap

mod_coefs <- mod1$coefficients
-mod_coefs[-1] * shap # PH coefficients

# Base rate
exp(mod_coefs[1])^(-shap) # or is this one?


# Other notes
# Doesnt converge with eps = 2, but will always with baseline model
survreg(Surv(t, eps == 2) ~ X1 + X2, data = dat, dist = "weibull")

