### Data generation for simulations ~ 


# Function to generate data:
# - two vars X and Z from multivar normal
# - Complete on Z, missing either MAR or MCAR or X
# - No censoring considered.


dat_gener <- function(N, # Sample size
                      X_type, # "contin" or "binary"
                      mus, # means of X and Z respectively
                      covmat, # Covariance matrix of X and Z
                      pars_logist, # Vector of intercept and slope for log reg (if X_type == binary)
                      mu_var, # Vector of mean and var for Z (if X_type == binary)
                      mech, # Specify "MCAR" or "MAR" 
                      p, # probability missingness; for MAR this is p when Z < 0; so 0.1 in MRR example
                      cause2, # Distribution for cause 2 times - either"weib" or "unif"
                      vals_t1, # vector of (a1, h1_0, B1_X, B1_Z)
                      vals_t2) { # vector of (a2, h2_0, B2_X, B2_Z)
  
  
  # Generate covariates
  if (X_type == "contin") {
    dat_covars <- data.frame(mvrnorm(n = N, mu = mus, Sigma = covmat)) %>%
      rename("X" = "X1", "Z" = "X2")
  } else {
    Z <- rnorm(n = N, mean = mu_var[1], sd = sqrt(mu_var[2]))
    pr <- 1 / (1 + exp(-(pars_logist[1] + pars_logist[2] * Z)))
    dat_covars <- data.frame(X = rbinom(n = N, 1, pr), Z)
  }
  
  
  # Read-in paramaters - Cause 1
  a1 <- vals_t1[1]
  h1_0 <- vals_t1[2]
  B1_X <- vals_t1[3]
  B1_Z <- vals_t1[4]
  
  lam1 <- with(dat_covars, h1_0 * exp(-(B1_X * X + B1_Z * Z) / a1))
  
  # Read-in paramaters - Cause 1
  a2 <- vals_t2[1]
  h2_0 <- vals_t2[2]
  B2_X <- vals_t2[3]
  B2_Z <- vals_t2[4]
  
  lam2 <- with(dat_covars, h2_0 * exp(-(B2_X * X + B2_Z * Z) / a2))
  
  
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
    
    dat <- ampute(dat_covars, prop = p, patterns = c(0, 1), mech = "MCAR")$amp
  } else {
    
    dat <- ampute(dat_covars, prop = p, patterns = c(0, 1), mech = "MAR")$amp
  }
  
  # Bring dataset together 
  final_dat <- dat %>%
    mutate(t = t, eps = as.factor(eps)) %>%
    
    # Append original (unimputed) covariate
    mutate(X_orig = dat_covars$X,
           Z_orig = dat_covars$Z) %>%
    
    # Compute cumulative hazards
    mutate(ev1 = as.numeric(eps == 1),
           ev2 = as.numeric(eps == 2)) %>%
    arrange(t) %>%
    mutate(H1 = nelsonaalen(., t, ev1),
           H2 = nelsonaalen(., t, ev2)) %>%
    
    # Add interaction terms
    mutate(H1_Z = H1 * Z,
           H2_Z = H2 * Z)
  
  return(final_dat)
}


# Two example of use
#dat <- dat_gener(N = 200, 
#                 X_type = "contin",
#                 mus = c(0, 0), 
#                 covmat = matrix(c(1, 0.25, 
#                                   0.25, 1), nrow = 2), 
#                 mech = "MAR", 
#                 p = 0.5, 
#                 cause2 = "weib", 
#                 vals_t1 = c(0.3, 1, 0.5, 0.5), 
#                 vals_t2 = c(1.7, 0.5, -0.5, 0.5)) 


#dat <- dat_gener(N = 200, 
#                 X_type = "binary", 
#                 pars_logist = c(0, -2), 
#                 mu_var = c(0, 1),  
#                 mech = "MAR", 
#                 p = 0.5, 
#                 cause2 = "weib", 
#                 vals_t1 = c(0.3, 1, 0.5, 0.5), 
#                 vals_t2 = c(1.7, 0.5, -0.5, 0.5)) 

print("hello world")

print("Where is this going")

print("Hein was here")

print('Liesbeth was here')
# Merge now 
