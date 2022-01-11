generate_dat_threecomp <- function(n,
                                   X_type,
                                   r,
                                   ev1_pars,
                                   ev2_pars,
                                   ev3_pars,
                                   rate_cens,
                                   mech = NULL,
                                   eta1 = NULL,
                                   p = NULL) {
  
  # Generate covariates
  covars <- gen_covars(n = n, X_type = X_type, r = r)
  
  # Generate the latent times
  T1 <- rweibull_KM(
    n = n, 
    alph = ev1_pars$a1, 
    lam = ev1_pars$h1_0 * exp(drop(c(ev1_pars$b1, ev1_pars$gamm1) %*% t(as.matrix(covars))))
  )
  T2 <- rweibull_KM(
    n = n, 
    alph = ev2_pars$a2, 
    lam = ev2_pars$h2_0 * exp(drop(c(ev2_pars$b2, ev2_pars$gamm2) %*% t(as.matrix(covars))))
  )
  T3 <- rweibull_KM(
    n = n, 
    alph = ev3_pars$a3, 
    lam = ev3_pars$h3_0 * exp(drop(c(ev3_pars$b3, ev3_pars$gamm3) %*% t(as.matrix(covars))))
  )
  
  t_tilde <- pmin(T1, T2, T3)
  eps <- data.table::fcase(
    t_tilde == T1, 1,
    t_tilde == T2, 2,
    t_tilde == T3, 3
  )
  admin_cens <- 10
  time <- pmin(t_tilde, admin_cens)
  delta <- as.numeric(time < admin_cens) * eps
  
  # Add regular censoring
  cens <- stats::rexp(n = n, rate = rate_cens)
  t <- pmin(time, cens)
  eps <- as.numeric(t < cens) * delta
  
  dat <- cbind.data.frame(t, eps, covars)
  dat <- dat[order(dat$t), ]
  
  # Add cumulative hazards and interactions
  dat$ev1 <- with(dat, as.numeric(eps == 1))
  dat$ev2 <- with(dat, as.numeric(eps == 2))
  dat$ev3 <- with(dat, as.numeric(eps == 3))
  
  dat$H1 <- nelsaalen_timefixed(dat, "t", "ev1")
  dat$H2 <- nelsaalen_timefixed(dat, "t", "ev2")
  dat$H3 <- nelsaalen_timefixed(dat, "t", "ev3")
  
  dat$H1_Z <- with(dat, H1 * Z)
  dat$H2_Z <- with(dat, H2 * Z)
  dat$H3_Z <- with(dat, H3 * Z)
  
  # Induce missings
  dat <- induce_missings(n, dat, p, mech, eta1) %>% 
    
    # Append original (unimputed) covariate
    dplyr::mutate(X_orig = X, X = X_miss) %>% 
    dplyr::select(-X_miss) %>%
    dplyr::mutate(eps = as.factor(eps))
  
  if (X_type == "binary") {
    dat <- dplyr::mutate(.data = dat, X = as.factor(X), X_orig = as.factor(X_orig))
  }
    
  return(dat)
}

setup_mstate_threecomp <- function(dat) {
  
  tmat <- mstate::trans.comprisk(K = 3)
  covs <- c("X", "Z")
  
  msdat <- mstate::msprep(
    time = c(NA, rep("t", 3)),
    status = with(dat, cbind(NA, eps == 1, eps == 2, eps == 3)),
    trans = tmat,
    data = dat,
    keep = covs
  )
  msdat_exp <- mstate::expand.covs(msdat, covs = covs, longnames = FALSE)
  ms_form <- reformulate(
    response = "Surv(time, status)",
    termlabels = c(paste0("X.", seq_len(3)), paste0("Z.", seq_len(3)), "strata(trans)")
  )
  cox_long <- survival::coxph(ms_form, data = msdat_exp)
  return(cox_long)
}
