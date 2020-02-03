###########################################
## Data generation as in MRR draft paper ## 
###########################################


# Generates following data for two competing risks:

# - Covariates X1, X2,..,XN from multivar normal
# - X1 can be binary, depending on X2 through logreg
# - Complete on X2,..,XN, missing either MAR or MCAR for X1
# - Missingness induced with mechanism from MRR draft paper
# - No censoring considered
# - T drawn from independent Weibulls
# - For now missingness mechanisms depend only on X2

# - PACKAGES NEEDED: MASS, tidyverse


dat_gener_MRR <- function(N, # Sample size
                          X_type, # "contin" or "binary"
                          mus, # means of MVN of covars
                          covmat, # Covariance matrix of X
                          pars_logist, # Vector of intercept and slope for log reg (if X_type == binary)
                          mech, # Specify "MCAR" or "MAR" 
                          p, # vector P miss if c(X2 > 0, X2 < 0); in MRR draft it was c(0.9, 0.1)
                          cause2, # Distribution for cause 2 times - either"weib" or "unif"
                          vals_t1, # vector of (a1, h1_0, B1_X, B1_X2)
                          vals_t2) { # vector of (a2, h2_0, B2_X, B2_X2)
                             
  
  # Generate covariates
  if (X_type == "contin") {
    dat_covars <- data.frame(mvrnorm(n = N, mu = mus, Sigma = covmat)) %>%
      rename_all(~ paste("X", 1:length(mus), sep = ""))
    
  } else {
    # For binary case, X1 depends on X2 via logreg
    dat_covars <- data.frame(mvrnorm(n = N, mu = mus, Sigma = covmat)) %>%
      rename_all(~ paste("X", 2:(length(mus) + 1), sep = ""))
    
    pr <- with(dat_covars, plogis(pars_logist[1] + pars_logist[2] * X2))
    dat_covars <- dat_covars %>%
      mutate(X1 = rbinom(n = N, 1, pr))
  }
  
  
  # Read-in paramaters - Cause 1
  a1 <- vals_t1[1]
  h1_0 <- vals_t1[2]
  B1_X1 <- vals_t1[3]
  B1_X2 <- vals_t1[4]
  
  lam1 <- with(dat_covars, h1_0 * exp(-(B1_X1 * X1 + B1_X2 * X2) / a1))
  
  # Read-in paramaters - Cause 1
  a2 <- vals_t2[1]
  h2_0 <- vals_t2[2]
  B2_X1 <- vals_t2[3]
  B2_X2 <- vals_t2[4]
  
  lam2 <- with(dat_covars, h2_0 * exp(-(B2_X1 * X1 + B2_X2 * X2) / a2))
  
  
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
    
    # This is MRR's MAR
    tab_X2 <- table(dat_covars$X2 > 0) 
    
    dat <- dat_covars %>% 
      mutate(miss_ind = ifelse(X2 > 0, rbinom(tab_X2[2], 1, p[1]), 
                               rbinom(tab_X2[1], 1, p[2])),
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
