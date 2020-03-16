
library(survival)
library(smcfcs)



# Example error -----------------------------------------------------------



set.seed(1)

# Sample size
n <- 10000

# Covariates
X <- rnorm(n, mean = 0, sd = 1)
Z <- rnorm(n, mean = 1, sd = 2)

# Induce approx. 30% MCAR missingness in X
miss_ind <- rbinom(n, 1, prob = .3)
X_miss <- ifelse(miss_ind == 1, NA, X)

# Generate times and bind data, no censoring
dat <- cbind.data.frame("time" = rexp(n, rate = exp(X + 0.5 * Z)),
                        "eps" = 1,
                        "X" = X_miss,
                        "Z" = Z)

# Sort by time
dat <- dat[order(dat$time), ]

# Impute X using smcfcs - will yield error 
imps <- smcfcs(dat, smtype = "coxph",
               smformula = "Surv(time, eps) ~ X + Z",
               method = c("", "", "norm", ""),
               m = 5, numit = 5)



# Why it gives the error --------------------------------------------------


# Run null model 
mod <- coxph(Surv(time, eps) ~ 1, data = dat)
H0 <- basehaz(mod, centered = F)

nrow(H0) # 9998, two time points are missing

# Rows ... and ... are missing
miss_time <- which(!(dat$time %in% H0$time))
dat[miss_time, ]


diff(dat[miss_time[1] + seq(-3, 3), "time"])

diff(tail(sort(dat[dat$time <= .11, ]$time)))

View(dat[order(dat$time), ])
