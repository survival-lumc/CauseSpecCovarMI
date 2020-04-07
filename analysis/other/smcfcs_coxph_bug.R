##*************************************##
## smcfcs impacted by default handling ##
##  of (near) ties in Survival package ##
##                                     ##
##    27/3/2020, Ed Bonneville         ##                            
##*************************************##


# Using R 3.6.2
library(survival) # version 3.1.11
library(smcfcs) # version 1.4.0


# Example error with simulated data ---------------------------------------


set.seed(1)

# Sample size
n <- 10000

# Covariates
X <- rnorm(n, mean = 0, sd = 1)
Z <- rnorm(n, mean = 1, sd = 2)

# Induce approx. 30% MCAR missingness in X
miss_ind <- rbinom(n, 1, prob = .3)
X_miss <- ifelse(miss_ind == 1, NA, X)

# Generate survival times and bind data, no censoring
dat <- cbind.data.frame("time" = rexp(n, rate = exp(X + 0.5 * Z)),
                        "eps" = 1,
                        "X" = X_miss,
                        "Z" = Z)

# Sort by time
dat <- dat[order(dat$time), ]

head(dat)

# Impute X using smcfcs - will yield error 
# not anymore from smcfcs 1.4.1!
imps <- smcfcs(dat, 
               smtype = "coxph",
               smformula = "Surv(time, eps) ~ X + Z",
               method = c("", "", "norm", ""),
               m = 5, 
               numit = 5)


# Why it gives the error --------------------------------------------------


# Run null model, get baseline hazard
mod <- coxph(Surv(time, eps) ~ 1, data = dat)
H0 <- basehaz(mod, centered = F)

nrow(H0) # 9996, four time points are missing

# Which ones are missing?
miss_time <- which(!(dat$time %in% H0$time))
dat[miss_time, ]

# Let's take a look at patient 117, and the patients after/before 
ind_pat117 <- which(rownames(dat) == 117)
dat[ind_pat117 + c(-1, 0, 1), ]

# The time step is small, order of 1e-8
diff(dat[ind_pat117 + c(-1, 0, 1), ]$time)

# The aeqSurv routine from survival will then consider this a tie,
# and the tied time-point will be equal to the earliest of the two


# The solution ------------------------------------------------------------


# Change coxph.control, where default is timefix = TRUE,
mod_fixed <- coxph(Surv(time, eps) ~ 1, data = dat,
                   control = coxph.control(timefix = FALSE))

H0_fixed <- basehaz(mod_fixed, centered = F)

nrow(H0_fixed) # everything is there!

