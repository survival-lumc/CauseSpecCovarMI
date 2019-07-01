############################################
## Data generation using unobserved var W ## 
############################################


# Generates following data for two competing risks:

# - Covariates X1, X2,..,XN from multivar normal
# - X1 can be binary, depending on X2 through logreg
# - Complete on X2,..,XN, missing either MAR or MCAR for X1
# - Missingness induced with logreg including W and standardised log(t)
# - No censoring considered
# - T drawn from independent Weibulls
# - For now missingness mechanisms depend only on X2

# - PACKAGES NEEDED: MASS, mice, tidyverse


dat_gener_W <- function(N, # Sample size
                        X_type, # "contin" or "binary"
                        mus, # means of MVN of covars
                        covmat, # Covariance matrix of X
                        pars_logist, # Vector of intercept and slope for log reg (if X_type == binary)
                        mech, # Specify "MCAR" or "MAR" 
                        pars_MAR, # vector (b_X2, b_t) for MAR mechanism
                        p, # probability missingness; for MAR this is p when Z < 0; so 0.1 in MRR example
                        cause2, # Distribution for cause 2 times - either"weib" or "unif"
                        vals_t1, # vector of (a1, h1_0, B1_X, B1_X2)
                        vals_t2) { # vector of (a2, h2_0, B2_X, B2_X2)
  
  
  # Generate covariates
  if (X_type == "contin") {
    dat_covars <- data.frame(mvrnorm(n = N, mu = mus, Sigma = covmat)) %>%
      rename_all(~ paste("X", 1:length(mus), sep = ""))

  } else {
    dat_covars <- data.frame(mvrnorm(n = N, mu = mus, Sigma = covmat)) %>%
      rename_all(~ paste("X", 2:(length(mus) + 1), sep = ""))
    
    # For binary case, X1 depends on X2 via logreg
    pr <- with(dat_covars, plogis(pars_logist[1] + pars_logist[2] * X2))
    dat_covars <- dat_covars %>%
      mutate(X1 = rbinom(n = N, 1, pr))
  }
  
  # Generate unobserved covariate
  W <- rnorm(N)
  
  # Read-in paramaters - Cause 1
  a1 <- vals_t1[1]
  h1_0 <- vals_t1[2]
  B1_X1 <- vals_t1[3]
  B1_X2 <- vals_t1[4]
  
  lam1 <- with(dat_covars, h1_0 * exp(-(B1_X1 * X1 + B1_X2 * X2 + W) / a1))
  
  # Read-in paramaters - Cause 1
  a2 <- vals_t2[1]
  h2_0 <- vals_t2[2]
  B2_X1 <- vals_t2[3]
  B2_X2 <- vals_t2[4]
  
  lam2 <- with(dat_covars, h2_0 * exp(-(B2_X1 * X1 + B2_X2 * X2 + W) / a2))
  
  
  # Generating event times with two independent weibulls:
  
  # For cause 1
  t1 <- rweibull(n = N, shape = a1, scale = lam1)
  
  # For cause 2
  if (cause2 == "unif") {
    t2 <- runif(n = N, min = 0, max = 0.5)
  } else {
    t2 <- rweibull(n = N, shape = a2, scale = lam2)
  }
  
  # Format times and event indicator
  t <- pmin(t1, t2)
  eps <- ifelse(t1 < t2, 1, 2)
  
  
  # Induce missingness
  if (mech == "MCAR") {
    
    dat <- dat_covars %>%
      mutate(miss_ind = rbinom(N, 1, p),
             X1 = ifelse(miss_ind == 1, NA, X1))  %>%
      select(X1, X2:paste("X", ncol(dat_covars), sep = "")) # Reorder 
    
  } else {
    
    # Compute log(t) and standardise
    t_stand <- (log(t) - mean(log(t))) / sd(log(t))
    
    # Define function to solve
    MAR_mech <- function(alph, beta, gam, p) {
      with(dat_covars, expr = {
        pr <- plogis(alph + beta * X2 + gam * t_stand) 
        mean(pr) - p
      })
    }
    
    # Get root
    alph <- uniroot(MAR_mech, interval = c(-10, 10), extendInt = "yes", 
                         beta = pars_MAR[1], gam = pars_MAR[2], p = p)$`root`
    
    # Induce missingness
    pr <- with(dat_covars, plogis(alph + pars_MAR[1] * X2 + pars_MAR[2] * t_stand))
    
    dat <- dat_covars %>% 
      mutate(miss_ind = rbinom(N, 1, pr),
             X1 = ifelse(miss_ind == 1, NA, X1)) %>%
      select(X1, X2:paste("X", ncol(dat_covars), sep = "")) # Reorder 
  }
  
  # Bring dataset together 
  final_dat <- dat %>%
    mutate(t = t, 
           eps = as.factor(eps),
           
           # Append original (unimputed) covariate
           X1_orig = dat_covars$X1) %>%
    
    # Compute cumulative hazards
    mutate(ev1 = as.numeric(eps == 1),
           ev2 = as.numeric(eps == 2)) %>%
    arrange(t) %>%
    mutate(H1 = nelsonaalen(., t, ev1),
           H2 = nelsonaalen(., t, ev2),
           
           # Interaction terms
           H1_X2 = H1 * X2,
           H2_X2 = H2 * X2) 
  
  return(final_dat)
}
