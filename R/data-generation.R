##***************************************##
## Support functions for generating data ##
##***************************************##


# Main --------------------------------------------------------------------


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
#' @inheritParams gen_cmprsk_times
#' 
#' @return Data-frame with missings induced
#' 
#' @export
generate_dat <- function(n,
                         X_type,
                         r,
                         ev1_pars,
                         ev2_pars,
                         rate_cens,
                         mech = NULL,
                         eta1 = NULL,
                         p = NULL,
                         mod_type = "latent") {
  
  # For checks
  X <- X_miss <- eps <- ev1 <- ev2 <- H1 <- H2 <- Z <- X_orig <- . <- NULL

  # Generate covariates
  dat_covars <- gen_covars(
    n = n,
    X_type = X_type,
    r = r
  )
  
  # Generate event times
  event_times <- gen_cmprsk_times(
    n,
    dat_covars, 
    ev1_pars,
    ev2_pars,
    rate_cens, 
    mod_type = mod_type
  )
  
  # Add event times
  dat_times <- dat_covars %>% 
    dplyr::mutate(t = event_times$t, eps = as.factor(event_times$eps))
  
  # Small check to provide if eta1 (when not MCAR)
  # This will also allow complete datasets with mech = NULL and eta1 = NULL
  if(!(is.null(mech)) && (mech != "MCAR") && is.null(eta1))
    stop("If mechanism is other than MCAR, you must provide an eta1 value.")
  
  # Induce missings and format
  dat <- induce_missings(n, dat_times, p, mech, eta1) %>% 
    
    # Append original (unimputed) covariate
    dplyr::mutate(X_orig = X, X = X_miss) %>% 
    dplyr::select(-X_miss) %>% # remove redundant variable
    
    # Compute individual event indicators
    dplyr::mutate(
      ev1 = as.numeric(eps == 1),
      ev2 = as.numeric(eps == 2)
    ) %>% 
    
    # Compute (marginal) cumulative hazards
    dplyr::mutate(
      H1 = nelsaalen_timefixed(., t, ev1),
      H2 = nelsaalen_timefixed(., t, ev2)
    ) %>% 
    
    # Compute interaction terms
    dplyr::mutate(
      H1_Z = H1 * Z,
      H2_Z = H2 * Z
    ) %>%
    dplyr::arrange(t) 

  # Convert to factors if binary
  if (X_type == "binary") {
    dat <- dat %>% 
      dplyr::mutate(
        X = as.factor(X),
        X_orig = as.factor(X_orig)
      )
  }
  
  return(dat)
}

gen_covars <- function(n, X_type, r) {
  
  # Compute true R needed if binary
  if (X_type == "binary") {
    r <- pbiserial_to_pearson(p = 0.5, r_pb = r)
    covmat <- matrix(c(1, r$r, r$r, 1), nrow = 2)
    
    dat_covars <- data.frame(MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = covmat)) %>%
      dplyr::rename_all(~ c("X", "Z")) %>%
      dplyr::mutate(X = ifelse(X <= r$x_cut, 0, 1))
    
  } else if (X_type == "continuous") {
    
    # MVN
    covmat <- matrix(c(1, r, r, 1), nrow = 2)
    
    dat_covars <- data.frame(MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = covmat)) %>%
      dplyr::rename_all(~ c("X", "Z"))
    
  } else if (X_type == "ordcat") { # Not relevant anymore
    
    # Ordered categorical
    Z <- stats::rnorm(n = n, mean = 0, sd = 1)
    
    # Baseline probabilities
    base_probs <- c(0.5, 0.25, 0.25)
    break_points <- c(-Inf, stats::qlogis(cumsum(base_probs)))
    
    # Add shift, covar effect of 1
    X_latent <- 1 * Z + stats::rlogis(n = n, location = 0, scale = 1)
    X <- cut(X_latent, breaks = break_points, ordered_result = T)
    dat_covars <- cbind.data.frame(X, X_latent, Z)
    
  } else stop("X_type should be either 'binary', 'continuous' or 'ordcat'")
  
}


#' @title Obtain necessary pearson to generate binary variable.
#' 
#' @param p Proportion of 1s in binary variable
#' @param r_pb Desired point-biserial correlation
#' 
#' @return Pearson r to use when generating MVN data.
#' 
#' @noRd
pbiserial_to_pearson <- function(p, r_pb) {
  
  x0 <- stats::qnorm(p, lower.tail = F) # cutoff
  h <- stats::dnorm(x0) # 'ordinate' = density/height of curve
  r <- r_pb * sqrt(p * (1 - p)) / h
  return(list("r" = r, "x_cut" = x0))
}


#' @title Generate competing risks times + event indictor (2 events)
#'
#' @param n Sample size.
#' @param dat Data frame containing n rows with X and Z as columns
#' @param ev1_pars Named list of ("a1", "h1_0", "b1", "gamm1")
#' @param ev2_pars Named list of ("a2", "h2_0", "b2", "gamm2")
#' @param rate_cens Rate of exponential distributiion for censoring.
#' If 0 means no censoring is applied.
#' @param mod_type Either "latent", weibull times generated from 
#' separate weibull distribution, or "total", times generated from 
#' sum of cause-specific hazards (using inverse transform method). For educational
#' purposes - both methods yield virtually same results.
#' 
#' @return Data-frame with missings induced
#' 
#' @export
gen_cmprsk_times <- function(n,
                             dat,
                             ev1_pars,
                             ev2_pars,
                             rate_cens,
                             mod_type = "latent") { 
  
  # Get model matrix 
  options(contrasts = rep("contr.treatment", 2)) 
  mod_mat <- stats::model.matrix(~ X + Z, dat = dat)
  
  # Compute rates
  lam1 <- ev1_pars$h1_0 * exp(mod_mat %*% c(0, ev1_pars$b1, ev1_pars$gamm1))
  lam2 <- ev2_pars$h2_0 * exp(mod_mat %*% c(0, ev2_pars$b2, ev2_pars$gamm2))

  if (mod_type == "latent") {
    
    t1 <- rweibull_KM(n = n, alph = ev1_pars$a1, lam = lam1)
    t2 <- rweibull_KM(n = n, alph = ev2_pars$a2, lam = lam2)
    
    # Take minimum of two, and compute indicator
    t <- pmin(t1, t2)
    eps <- ifelse(t1 < t2, 1, 2)
    
  } else if (mod_type == "total") {
    
    # Use inverse transform sampling
    t <- invtrans_weib(
      n = n, 
      alph1 = ev1_pars$a1, 
      lam1 = lam1, 
      alph2 = ev2_pars$a2, 
      lam2 = lam2
    )
    
    # Determine which event occured
    haz_ev1 <- haz_weib(alph = ev1_pars$a1, lam = lam1, t = t)
    haz_ev2 <- haz_weib(alph = ev2_pars$a2, lam = lam2, t = t)
    event <- stats::rbinom(n = n, size = 1, prob = haz_ev1 / (haz_ev1 + haz_ev2))
    eps <- ifelse(event == 1, 1, 2)
  }
  
  # Add censoring
  if (0 < rate_cens) {
    cens <- stats::rexp(n = n, rate = rate_cens)
    eps <- ifelse(cens < t, 0, eps)
    t <- pmin(cens, t)
    
    # Also artificially censor at 10 years
    eps <- ifelse(t >= 10, 0, eps)
    t <- pmin(t, 10)
  }
  
  return(cbind.data.frame(t, eps))
}




invtrans_weib <- function(n, alph1, lam1, alph2, lam2) {
  
  t_tilde <- NULL
  
  # Define cdf - U first
  cdf_U <- function(t, alph1, lam1, alph2, lam2, U) {
    F_min_U <- 1 - exp(-(lam1 * t^alph1 + lam2 * t^alph2)) - U
    return(F_min_U)
  }
  
  # Generate u and store in df
  u <- stats::runif(n)
  dat_roots <- cbind.data.frame(u, alph1, lam1, alph2, lam2)
  
  # Generate uniform values, and find roots - better to use later rstpm2::vuniroot()
  samps <- dat_roots %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      t_tilde = stats::uniroot(
        cdf_U,
        interval = c(.Machine$double.eps,
                     1000),
        extendInt = "yes",
        U = u,
        alph1 = alph1,
        lam1 = lam1,
        alph2 = alph2,
        lam2 = lam2
      )$`root`
    ) %>%
    dplyr::pull(t_tilde)
  
  return(samps)
}


#' @title Sample from Weibull distirbution in K&M parametrisation.
#' 
#' @param alph Shape parameter of weibull distribution.
#' @param lam Rate parameter of lambda distirbution.
#' @param n sample size
#' 
#' @return n samples from weibull distribution.
rweibull_KM <- function(n, alph, lam) {
  
  samp <- (-log(1 - stats::runif(n)) / lam)^(1 / alph)
  return(samp)
}


#' @title Induce missingess.
#' 
#' @param n Sample size.
#' @param dat Data frame containing X and Z
#' @param mech Missingness mechanism, one of:
#' "MAR_GEN", "MAR", "MNAR", "MCAR"
#' @param eta1 Only necessary for mech != "MCAR" - degree/direction
#' of assocation between variable responsible for missingness
#' and probability of missingness
#' @param p Proportion of missing values in X.
#' 
#' @return Dataset with missingness induced.
#' 
#' @noRd
induce_missings <- function(n, dat, p, mech, eta1) {
  
  X <- miss_ind <- NULL

  if (is.null(mech) | is.null(p)) {
    pr <- 0
 
  } else if (mech == "MCAR") {
    pr <- p
    
  } else if (mech == "MAR") {
    pr <- logreg_missings(p, eta1, covar = dat$Z)

  } else if (mech == "MNAR") {
    
    # Check if X ordered cat
    if (length(levels(dat$X)) > 2) {
      pr <- logreg_missings(p, eta1, covar = dat$X_latent)
    } else {
      pr <- logreg_missings(p, eta1, covar = dat$X)
    }
    
  } else if (mech == "MAR_GEN") {
    pr <- logreg_missings(p, eta1, covar = scale(log(dat$t)))
  }
  
  dat <- dat %>%
    dplyr::mutate(
      miss_ind = stats::rbinom(n, 1, pr),
      X_miss = ifelse(miss_ind == 1, NA, X)
    )
  
  return(dat)
}


#' @title Calculate probability of missingness using 
#' logistic regression.
#' 
#' @param p Proportion of missing values in X.
#' @param eta1 Only necessary for mech != "MCAR" - degree/direction
#' of assocation between variable responsible for missingness
#' and probability of missingness
#' @param covar Variable responsible for missingness, e.g. "Z"
#' 
#' @return Probability vector to generate missings.
#' 
#' @noRd
logreg_missings <- function(p, eta1, covar) {
  
  intercept_solve <- function(eta0, eta1, p) {
    pr <- stats::plogis(eta0 + eta1 * covar) 
    return(mean(pr) - p)  
  }
  
  eta0 <- stats::uniroot(
    intercept_solve, 
    interval = c(-25, 25), 
    extendInt = "yes", 
    eta1 = eta1, 
    p = p
  )$`root` 
  
  # Induce missingness
  pr <- stats::plogis(eta0 + eta1 * covar)
  
  return(pr)
}

#' Compute Nelson Aalen estimate for simulated data
#' 
#' Version of \code{mice::nelsonaalen} but with the
#' \code{control = survival::coxph.control(timefix = FALSE)} argument added
#' to the \code{survival::coxph} call. This is to deal with floating point errors
#' within the simulation procedure, where two or more simulated time points
#' with minimal difference between would normally be considered as a tie.
#' 
#' @param dat Data frame with time and status variables
#' @param timevar Column name of time variable
#' @param statusvar Column name of (numeric) status variable
#' @param timefix Logical, whether to apply timefix or not. If TRUE it corresponds 
#' to using \code{mice::nelsonaalen}
#' 
#' @export
nelsaalen_timefixed <- function(dat,
                                timevar,
                                statusvar,
                                timefix = FALSE) {
  
  timevar <- as.character(substitute(timevar))
  statusvar <- as.character(substitute(statusvar))
  time <- dat[, timevar]
  status <- dat[, statusvar]
  mod <- survival::coxph(
    Surv(time, status) ~ 1, 
    control = survival::coxph.control(timefix = timefix)
  )
  hazard <- survival::basehaz(mod)
  idx <- match(time, hazard[, "time"])
  return(hazard[idx, "hazard"])
}
