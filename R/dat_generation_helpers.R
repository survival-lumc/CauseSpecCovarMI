##***************************************##
## Support functions for generating data ##
##***************************************##


# Set-up ------------------------------------------------------------------


options(readr.num_columns = 0) # suppress read.csv message

# Add global variable (for package errors)
globalVariables(c(names(
  readr::read_csv(system.file("testdata",
                              "test_data.csv",
                              package = "SimsCauseSpecCovarMiss"))
), "object", "res", "times", "variance", "warn", ".", "X_miss"))


# Main --------------------------------------------------------------------


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
  #' @importFrom magrittr `%>%` 
  #' @importFrom rlang .data
  #' 
  #' @return Data-frame with missings induced
  #' 
  #' @export

  # Compute true R needed if binary
  if (X_type == "binary") {
    r <- pbiserial_to_pearson(p = 0.5, r_pb = r)

    covmat <- matrix(c(1, r$r,
                       r$r, 1), nrow = 2)

    dat_covars <- data.frame(
      MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = covmat)
    ) %>%
      dplyr::rename_all(~ c("X", "Z")) %>%
      dplyr::mutate(X = ifelse(X <= r$x_cut, 0, 1))
  } else {
    covmat <- matrix(c(1, r,
                       r, 1), nrow = 2)

    dat_covars <- data.frame(
      MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = covmat)
    ) %>%
      dplyr::rename_all(~ c("X", "Z"))
  }
  
  # Generate event times
  event_times <- gen_cmprsk_times(n,
                                  dat_covars, 
                                  ev1_pars,
                                  ev2_pars,
                                  rate_cens, 
                                  mod_type = mod_type)
  # Add event times
  dat_times <- dat_covars %>% 
    dplyr::mutate(t = event_times$t, 
                  eps = as.factor(event_times$eps))
  
  # Small check to provide if eta1 (when not MCAR)
  # This will also allow complete datasets with mech = NULL and eta1 = NULL
  if(!(is.null(mech)) && (mech != "MCAR") && is.null(eta1))
    stop("If mechanism is other than MCAR, you must provide an eta1 value.")
  
  # Induce missings and format
  dat <- induce_missings(n, dat_times,
                         p, mech, eta1) %>% 
    
    # Append original (unimputed) covariate
   dplyr::mutate(X_orig = X,
                 X = X_miss) %>% 
    dplyr::select(-X_miss) %>% # remove redundant variable
    
    # Compute individual event indicators
    dplyr::mutate(ev1 = as.numeric(eps == 1),
                  ev2 = as.numeric(eps == 2)) %>% 
    
    # Compute (marginal) cumulative hazards
    dplyr::mutate(H1 = nelsaalen_timefixed(., t, ev1),
                  H2 = nelsaalen_timefixed(., t, ev2)) %>% 
    
    # Compute interaction terms
    dplyr::mutate(H1_Z = H1 * Z,
                  H2_Z = H2 * Z) %>%
    dplyr::arrange(t) 

  # Convert to factors if binary
  if (X_type == "binary") {
    dat <- dat %>% 
      dplyr::mutate(X = as.factor(dat$X),
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
  
  x0 <- stats::qnorm(p, lower.tail = F) # cutoff
  h <- stats::dnorm(x0) # 'ordinate' = density/height of curve
  r <- r_pb * sqrt(p * (1 - p)) / h
  return(list("r" = r, "x_cut" = x0))
}


gen_cmprsk_times <- function(n,
                             dat,
                             ev1_pars,
                             ev2_pars,
                             rate_cens,
                             mod_type = "latent") { 
  
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
  #' sum of cause-specific hazards (using inverse transform method)
  #' 
  #' @return Data-frame with missings induced
  
  # Rates from cause specific hazards
  lam1 <- with(
    dat,
    ev1_pars$h1_0 * exp((ev1_pars$b1 * X + ev1_pars$gamm1 * Z))
  )
  
  lam2 <- with(
    dat,
    ev2_pars$h2_0 * exp((ev2_pars$b2 * X + ev2_pars$gamm2 * Z))
  )
  
  if (mod_type == "latent") {
    
    t1 <- rweibull_KM(n = n, alph = ev1_pars$a1, lam = lam1)
    t2 <- rweibull_KM(n = n, alph = ev2_pars$a2, lam = lam2)
    
    # Take minimum of two, and compute indicator
    t <- pmin(t1, t2)
    eps <- ifelse(t1 < t2, 1, 2)
    
  } else if (mod_type == "total") {
    
    # n = 1 is misleading: this is vectorised so we actually obtain n samples
    t <- invtrans_weib(
      n = 1, 
      alph1 = ev1_pars$a1, 
      lam1 = lam1, 
      alph2 = ev2_pars$a2, 
      lam2 = lam2
    )
    
    # Determine which event occured
    haz_ev1 <- haz_weib(ev1_pars$a1, lam1, t)
    haz_ev2 <- haz_weib(ev2_pars$a2, lam2, t)
    event <- stats::rbinom(n, 1, prob = haz_ev1 / (haz_ev1 + haz_ev2))
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




invtrans_weib <- Vectorize(function(n, alph1, lam1, alph2, lam2) {
  
  # Define cdf - U first
  cdf_U <- function(t, U) {
    F_min_U <- 1 - exp(-(lam1 * t^alph1 + lam2 * t^alph2)) - U
    return(F_min_U)
  }
  
  # Generate uniform values, and find roots
  samps <- sapply(stats::runif(n), function(u) {
    root <- stats::uniroot(
      cdf_U, interval = c(.Machine$double.eps, 1000), 
      extendInt = "yes", U = u
    )$`root`
    return(root)
  })
  
  return(samps)
})




rweibull_KM <- function(n, alph, lam) {
  
  #' @title Sample from weibull in K&M parametrisation.
  #' 
  #' @param alph Shape parameter of weibull distribution.
  #' @param lam Rate parameter of lambda distirbution.
  #' @param n sample size
  #' 
  #' @return n samples from weibull distribution.
  
  samp <- (-log(1 - stats::runif(n)) / lam)^(1 / alph)
  return(samp)
}


induce_missings <- function(n, dat, p, mech, eta1) {
  
  #' @title Induce missingess.
  #' 
  #' @param n Sample size.
  #' @param dat Data fram containing X and Z
  #' @param mech Missingness mechanism, one of:
  #' "MAR_GEN", "MAR", "MNAR", "MCAR"
  #' @param eta1 Only necessary for mech != "MCAR" - degree/direction
  #' of assocation between variable responsible for missingness
  #' and probability of missingness
  #' @param p Proportion of missing values in X.
  #' 
  #' @return Dataset with missingness induced.

  if (is.null(mech) | is.null(p)) {
    pr <- 0
 
  } else if (mech == "MCAR") {
    pr <- p
    
  } else if (mech == "MAR") {
    pr <- logreg_missings(p, eta1, covar = dat$Z)

  } else if (mech == "MNAR") {
    pr <- logreg_missings(p, eta1, covar = dat$X)
    
  } else if (mech == "MAR_GEN") {
    pr <- logreg_missings(p, eta1, covar = scale(log(dat$t)))
  }
  
  dat <- dat %>%
    dplyr::mutate(miss_ind = stats::rbinom(n, 1, pr),
                  X_miss = ifelse(miss_ind == 1, NA, X))
  
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
  #' @return Probability vector to generate missings.
  
  intercept_solve <- function(eta0, eta1, p) {
    pr <- stats::plogis(eta0 + eta1 * covar) 
    return(mean(pr) - p)  
  }
  
  eta0 <- stats::uniroot(intercept_solve, 
                         interval = c(-25, 25), 
                         extendInt = "yes", 
                         eta1 = eta1, p = p)$`root` 
  
  # Induce missingness
  pr <- stats::plogis(eta0 + eta1 * covar)
  
  return(pr)
}


nelsaalen_timefixed <- function(dat,
                                timevar,
                                statusvar,
                                timefix = FALSE) {
  
  #' @importFrom survival Surv strata
  
  timevar <- as.character(substitute(timevar))
  statusvar <- as.character(substitute(statusvar))
  time <- dat[, timevar]
  status <- dat[, statusvar]
  mod <- survival::coxph(Surv(time, status) ~ 1, 
                         control = survival::coxph.control(timefix = timefix))
  hazard <- survival::basehaz(mod)
  idx <- match(time, hazard[, "time"])
  return(hazard[idx, "hazard"])
}
