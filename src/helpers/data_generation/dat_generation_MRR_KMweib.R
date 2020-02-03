############################################
## Data generation as in MRR draft paper  ## 
##          with KM weib param            ## 
############################################


# Generates following data for two competing risks:

# - Covariates X1, X2,..,XN from multivar normal
# - X1 can be binary, depending on X2 through logreg
# - Complete on X2,..,XN, missing either MAR or MCAR for X1
# - Missingness induced with mechanism from MRR draft paper
# - No censoring considered
# - T drawn from independent Weibulls
# - For now missingness mechanisms depend only on X2

# - PACKAGES NEEDED: MASS, tidyverse

# Alternative parametrisation
my_rweibull <- function(n, alph, lam) {
  return((-log(1 - runif(n)) / lam)^(1 / alph))
}

haz_weib <- Vectorize(function(alph, lam, t) {
  return(alph * lam * t^(alph - 1))
})

dens_weib <- Vectorize(function(alph, lam, t) {
  haz_weib(alph, lam, t) * exp(-lam * t^alph)
})

cumhaz_weib <- Vectorize(function(alph, lam, t) {
  lam * t^alph
})

# Take input vector of
gen_surv_weib <- Vectorize(function(cumhaz1, cumhaz2) {
  exp(-(cumhaz1 + cumhaz2))
})


# Make a general cumulative incidence function
cuminc_weib <- function(alph_ev, lam_ev, alph_comp, lam_comp, t) {
  
  prod <-  function(t) {
    haz_weib(alph_ev, lam_ev, t) * gen_surv_weib(cumhaz_weib(alph_ev, lam_ev, t),
                                                 cumhaz_weib(alph_comp, lam_comp, t))
  }
  
  ci_func <- Vectorize(function(upp) {
    integrate(prod, lower = 0, upper = upp)$value
  })
  
  return(ci_func(t))
}


#ci1 <- cuminc_weib(alph_ev = 0.3, lam_ev = 1, 
#                   alph_comp = 1.7, lam_comp = 1.2,
#                   t = dat$t)

#ci2 <- cuminc_weib(alph_ev = 1.7, lam_ev = 1.2, 
#                   alph_comp = 0.3, lam_comp = 1,
#                   t = dat$t)




dat_gener_MRR_KM <- function(N, # Sample size
                             X_type, # "contin" or "binary"
                             mus, # means of MVN of covars
                             covmat, # Covariance matrix of X
                             pars_logist, # Vector of intercept and slope for log reg (if X_type == binary)
                             mech, # Specify "MCAR" or "MAR" 
                             p, # vector P miss if c(X2 > 0, X2 < 0); in MRR draft it was c(0.9, 0.1)
                             cause2, # Distribution for cause 2 times - either"weib" or "unif"
                             vals_t1, # vector of (a1, h1_0, B1_X, B1_X2)
                             vals_t2, # vector of (a2, h2_0, B2_X, B2_X2)
                             rate_cens) { 
                          
  
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
  
  lam1 <- with(dat_covars, h1_0 * exp((B1_X1 * X1 + B1_X2 * X2))) # what about minus?
  
  # Read-in paramaters - Cause 1
  a2 <- vals_t2[1]
  h2_0 <- vals_t2[2]
  B2_X1 <- vals_t2[3]
  B2_X2 <- vals_t2[4]
  
  lam2 <- with(dat_covars, h2_0 * exp((B2_X1 * X1 + B2_X2 * X2)))
  
  
  # Generating event times with two independent weibulls:
  
  # For cause 1
  t1 <- my_rweibull(n = N, alph = a1, lam = lam1)
  
  # For cause 2
  if (cause2 == "unif") {
    t2 <- runif(n = N, min = 0, max = 0.5)
  } else {
    t2 <- my_rweibull(n = N, alph = a2, lam = lam2)
  }
  
  # Format times and event indicator
  t <- pmin(t1, t2)
  eps <- ifelse(t1 < t2, 1, 2)
  
  # Censoring 
  if (rate_cens != 0) {
    cens <- rexp(n = N, rate = rate_cens)
    eps <- ifelse(cens < t, 0, eps)
  }
  
  
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
    
    # Compute hazards
    #mutate(haz1 = haz_weib(alph = a1, lam = lam1, t = t),
    #       haz2 = haz_weib(alph = a2, lam = lam2, t = t)) %>% 
    arrange(t) %>%
    mutate(H1 = nelsonaalen(., t, ev1),
           H2 = nelsonaalen(., t, ev2),
           
           # Interaction terms
           H1_X2 = H1 * X2,
           H2_X2 = H2 * X2) 
  return(final_dat)
}


# Comment out the bottom bit

# test run
