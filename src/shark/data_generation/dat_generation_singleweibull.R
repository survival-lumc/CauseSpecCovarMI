##########################################
## Data generation using single weibull ## 
##########################################


# Generates following data for two competing risks:

# - Covariates X1, X2,..,XN from multivar normal
# - X1 can be binary, depending on X2 through logreg
# - Complete on X2,..,XN, missing either MAR or MCAR for X1
# - Missingness induced with mechanism from MRR draft paper
# - No censoring considered
# - T drawn from single Weibull
# - For now missingness mechanisms depend only on X2

# - PACKAGES NEEDED: MASS, tidyverse


### Support functions ~


# 1) Generates n samples of Weibull(a, b) parametrisation in Klein & Moeschberger (p38 overview)
my_rweibull <- function(n, alph, lam) {
  return((-log(1 - runif(n)) / lam)^(1 / alph))
}


# 2) PDF Weibull in K&M parametrisation
pdf_weib <- function(t, alph, lam) {
  haz <- alph * lam * t^(alph - 1)
  sur <- exp(-lam * t^alph)
  return(haz * sur)
}


# 3) Hazard function of Weibull in K&M parametrisation
haz_weib <- function(alph, lam, t) {
  return(alph * lam * t^(alph - 1))
}


# 4) Inverse-transform using uniroot
invtrans_weib <- Vectorize(function(n, alph1, lam1, alph2, lam2) {
  
  # Define cdf - U first
  cdf_U <- function(t, U) {
    F_min_U <- 1 - exp(-(lam1 * t^alph1 + lam2 * t^alph2)) - U
    return(F_min_U)
  }
  
  # Generate uniform values, and find roots
  samps <- sapply(runif(n), function(u) {
    root <- uniroot(cdf_U, interval = c(.Machine$double.eps, 1000), 
                    extendInt = "yes", U = u)$`root`
    return(root)
  })
  
  return(samps)
})


### Function to generate data ~ 


# - Covariates X1, X2,..,XN from multivar normal
# - Complete on X2,..,XN, missing either MAR or MCAR or X1
# - No censoring considered.
# - For now missingness mechanisms depend only on X2


dat_gener_singlweib <- function(N, # Sample size
                                X_type, # "contin" or "binary"
                                mus, # means of MVN of covars
                                covmat, # Covariance matrix of X
                                pars_logist, # Vector of intercept and slope for log reg (if X_type == binary)
                                mech, # Specify "MCAR" or "MAR" 
                                p, # probability missingness; for MAR this is p when Z < 0; so 0.1 in MRR example
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
  
  
  # Generating event times with single weibull:
  
  # Change parametrisation first
  lam1 <- lam1^(-a1)
  lam2 <- lam2^(-a2)
  
  # n = 1 is misleading: this is vectorised so we actually obtain N samples
  t <- invtrans_weib(n = 1, alph1 = a1, lam1 = lam1, alph2 = a2, lam2 = lam2)
  
  # Determine which event occured
  ratio <- haz_weib(a1, lam1, t) / (haz_weib(a1, lam1, t) + haz_weib(a2, lam2, t))
  event <- rbinom(N, 1, ratio)
  eps <- ifelse(event == 1, 1, 2)
  
  
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
      mutate(miss_ind = ifelse(X2 > 0, rbinom(tab_X2[2], 1, p), 
                               rbinom(tab_X2[1], 1, 1 - p)),
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
