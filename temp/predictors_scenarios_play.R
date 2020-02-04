
dato <- cbind.data.frame(ric = as.numeric((dat$ric)) - 1,
                         age = as.numeric(as.character(dat_orig$age_allo1)))

hist(dato$age)

dat$dato <- dato[complete.cases(dato), ]


table(dato$ric)/nrow(dato)

mod <- glm(ric ~ age, family = binomial, data = dato)

eta <- cbind(1, dato$age) %*% mod$coefficients

y <- rbinom(n=nrow(dato), size=1, prob=plogis(eta))

table(y)/length(y)


karno <- as.numeric(dat_orig$KARNOFSK_allo1)
hist(karno)

library(MASS)

sigma <- cor(dato)

lol <- as.data.frame(mvrnorm(n = nrow(dato), mu = c(mean(dato$age), ),
        Sigma = sigma))

missings <- miss_var_summary(dat_orig)

# Scenarios
beta <- c(0, 0.5, 1)
samp_size <- c(500, 2000)
miss <- c("low", "high")
X1_type <- c("contin", "binary")
miss_scen <- c("MCAR", "MAR_MRR", "MNAR",
               "MAR_gen")

scens <- expand.grid("beta" = beta, 
                     "miss" = miss, 
                     "X1_type" = X1_type, 
                     "scen" = miss_scen,
                     "n" = samp_size)

scens
