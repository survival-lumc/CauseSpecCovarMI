##***************************************##
## Support functions for generating data ##
##***************************************##

# Change import from dplyr and tidyr to specific functions

generate_dat <- function(n,
                         X_type,
                         r,
                         ev1_pars,
                         ev2_pars,
                         rate_cens,
                         mech,
                         eta1,
                         p) {

  #' @title Generates data.
  #'
  #' @param n Sample size.
  #' @param X_type Either "binary", or "contin"
  #' @param r Desired correlation between X and Z.
  #' @param ev1_pars Named list of ("a1" =, "h1_0", "b1", "gamm1")
  #' @param ev2_pars Named list of ("a2", "h2_0", "b2", "gamm2")
  #' @param rate_cens Rate of exponential distributiion for censoring.
  #' If 0 means no censoring is applied.
  #' @param mech Missingness mechanism, one of:
  #' "MAR_GEN", "MAR", "MNAR", "MCAR"
  #' @param eta1 Only necessary for mech != "MCAR" - degree/direction
  #' of assocation between variable responsible for missingness
  #' and probability of missingness
  #' @param p Proportion of missing values in X.
  #' 
  #' @importFrom magrittr `%>%`
  #' @importFrom MASS mvrnorm
  #' @importFrom mice nelsonaalen
  #' @import dplyr
  #' @import tidyr
  #' @import stringr
  #' @import ggplot2
  #' 
  #' @return Data-frame with missings induced

  # Compute true R needed if binary
  if (X_type == "binary") {
    r <- pbiserial_to_pearson(p = 0.5, r_pb = r)

    covmat <- matrix(c(1, r$r,
                       r$r, 1), nrow = 2)

    dat_covars <- data.frame(
      mvrnorm(n = n, mu = c(0, 0), Sigma = covmat)
    ) %>%
      rename_all(~ c("X", "Z")) %>%
      mutate(X = ifelse(X <= r$x_cut, 0, 1))
  } else {
    covmat <- matrix(c(1, r,
                       r, 1), nrow = 2)

    dat_covars <- data.frame(
      mvrnorm(n = n, mu = c(0, 0), Sigma = covmat)
    ) %>%
      rename_all(~ c("X", "Z"))
  }
  
  # Generate event times
  event_times <- gen_cmprsk_times(n,
                                  dat_covars, 
                                  ev1_pars,
                                  ev2_pars,
                                  rate_cens)
  # Add event times
  dat_times <- dat_covars %>% 
    mutate(t = event_times$t, 
           eps = as.factor(event_times$eps))
  
  # Induce missings and format
  dat <- induce_missings(n, dat_times,
                         p, mech, eta1) %>% 
    
    # Append original (unimputed) covariate
    mutate(X_orig = X, # 
           X = X_miss) %>% 
    select(-X_miss) %>% #remove redundant variable
    
    # Compute cumulative hazards
    mutate(ev1 = as.numeric(eps == 1),
           ev2 = as.numeric(eps == 2)) %>% 
    mutate(H1 = nelsonaalen(., t, ev1),
           H2 = nelsonaalen(., t, ev2),
           
           # Interaction terms
           H1_Z = H1 * Z,
           H2_Z = H2 * Z) %>%
    
    # Compute hazards - dleet if not needed
    #mutate(haz1 = haz_weib(alph = a1, lam = lam1, t = t),
    #       haz2 = haz_weib(alph = a2, lam = lam2, t = t)) %>% 
    arrange(t) 
  
  # Convert to factors if binary
  if (X_type == "binary") {
    dat <- dat %>% 
      mutate(X = as.factor(dat$X),
             X_orig = as.factor(dat$X_orig))
  }
  
  return(dat)
}


pbiserial_to_pearson <- function(p, r_pb) {
  
  #' @title Obtain necessary pearson to generate binary variable.
  #' 
  #' @param p Proportion of 1s in binary variable
  #' @param r_pb Desired point-biserial correlation
  #' 
  #' @return Pearson r to use when generating MVN data.
  
  x0 <- qnorm(p, lower.tail = F) # cutoff
  h <- dnorm(x0) # 'ordinate' = density/height of curve
  r <- r_pb * sqrt(p * (1 - p)) / h
  return(list("r" = r, "x_cut" = x0))
}


gen_cmprsk_times <- function(n,
                             .data,
                             ev1_pars,
                             ev2_pars,
                             rate_cens) {
  
  #' @title Generate competing risks times + event indictor (2 events)
  #'
  #' @param n Sample size.
  #' @param .data Data frame containing n rows with X and Z as columns
  #' @param ev1_pars Named list of ("a1", "h1_0", "b1", "gamm1")
  #' @param ev2_pars Named list of ("a2", "h2_0", "b2", "gamm2")
  #' @param rate_cens Rate of exponential distributiion for censoring.
  #' If 0 means no censoring is applied.
  #' 
  #' @return Data-frame with missings induced
  
  # Generate time to events
  lam1 <- with(
    .data,
    ev1_pars$h1_0 * exp((ev1_pars$b1 * X + ev1_pars$gamm1 * Z))
  )
  
  t1 <- rweibull_KM(n = n, alph = ev1_pars$a1, lam = lam1)
  
  lam2 <- with(
    .data,
    ev2_pars$h2_0 * exp((ev2_pars$b2 * X + ev2_pars$gamm2 * Z))
  )
  
  t2 <- rweibull_KM(n = n, alph = ev2_pars$a2, lam = lam2)
  
  # Take minimum of two, and compute indicator
  t <- pmin(t1, t2)
  eps <- ifelse(t1 < t2, 1, 2)
  
  # Add censoring
  if (0 < rate_cens) {
    cens <- rexp(n = n, rate = rate_cens)
    eps <- ifelse(cens < t, 0, eps)
    t <- pmin(cens, t)
  }
  
  return(cbind.data.frame(t, eps))
}


rweibull_KM <- function(n, alph, lam) {
  
  #' @title Sample from weibull in K&M parametrisation.
  #' 
  #' @param alph Shape parameter of weibull distribution.
  #' @param lam Rate parameter of lambda distirbution.
  #' 
  #' @return n samples from weibull distribution.
  
  samp <- (-log(1 - runif(n)) / lam)^(1 / alph)
  return(samp)
}


induce_missings <- function(n, .data, p, mech, eta1) {
  
  #' @title Induce missingess.
  #' 
  #' @param n Sample size.
  #' @param .data Data fram containing X and Z
  #' @param mech Missingness mechanism, one of:
  #' "MAR_GEN", "MAR", "MNAR", "MCAR"
  #' @param eta1 Only necessary for mech != "MCAR" - degree/direction
  #' of assocation between variable responsible for missingness
  #' and probability of missingness
  #' @param p Proportion of missing values in X.
  #' 
  #' @return 
  
  if (mech == "MCAR") {
    dat <- .data %>%
      mutate(miss_ind = rbinom(n, 1, p),
             X_miss = ifelse(miss_ind == 1, NA, X))
    
  } else if (mech == "MAR") {
    pr <- logreg_missings(p, eta1, covar = .data$Z)
    
    dat <- .data %>%
      mutate(miss_ind = rbinom(n, 1, pr),
             X_miss = ifelse(miss_ind == 1, NA, X))
    
  } else if (mech == "MNAR") {
    pr <- logreg_missings(p, eta1, covar = .data$X)
    
    dat <- .data %>%
      mutate(miss_ind = rbinom(n, 1, pr),
             X_miss = ifelse(miss_ind == 1, NA, X))
    
  } else if (mech == "MAR_GEN") {
    pr <- logreg_missings(p, eta1, covar = scale(log(.data$t)))
    
    dat <- .data %>%
      mutate(miss_ind = rbinom(n, 1, pr),
             X_miss = ifelse(miss_ind == 1, NA, X))
    
  }
  
  return(dat)
}


logreg_missings <- function(p, eta1, covar) {
  
  #' @title Calculate probability of missingness using 
  #' logistic regression.
  #' 
  #' @param p Proportion of missing values in X.
  #' @param eta1 Only necessary for mech != "MCAR" - degree/direction
  #' of assocation between variable responsible for missingness
  #' and probability of missingness
  #' @param covar Variable responsible for missingness, e.g. "Z"
  #' 
  #' @return 
  
  intercept_solve <- function(eta0, eta1, p) {
    pr <- plogis(eta0 + eta1 * covar) 
    return(mean(pr) - p)  
  }
  
  eta0 <- uniroot(intercept_solve, 
                  interval = c(-25, 25), 
                  extendInt = "yes", 
                  eta1 = eta1, p = p)$`root` # change eta1 to strength
  
  # Induce missingness
  pr <- plogis(eta0 + eta1 * covar)
  
  return(pr)
}

