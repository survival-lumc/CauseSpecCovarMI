##************************##
## Testing mstate helpers ##
##************************##


setup_mstate <- function(dat) {
  
  #' @title Set-up mstate model (pre-predicting)
  #' 
  #' @description ...
  #' 
  #' @inheritParams get_predictor_mats
  #' 
  #' @importFrom mstate msprep trans.comprisk expand.covs msfit probtrans
  #' @importFrom survival coxph Surv strata survfit
  #' 
  #' @return Long cox model
  #' 
  #' @export
  
  # Set up transition matrix 
  tmat <- mstate::trans.comprisk(2, c("Rel", "NRM"))
  covs <- c("X", "Z")
  
  # Long format
  dat_msprepped <- mstate::msprep(time = c(NA, "t", "t"),
                                  status = c(NA, "ev1", "ev2"), 
                                  data = dat,
                                  trans = tmat,
                                  keep = covs) 
                          
  # Expand covariates
  dat_expanded <- mstate::expand.covs(dat_msprepped, covs,
                                      append = TRUE, longnames = F)
  
  # Fit long cox model (both transitions)
  cox_long <- survival::coxph(Surv(time, status) ~ 
                                X.1 + Z.1 + # Trans == 1
                                X.2 + Z.2 + # Trans == 2
                                strata(trans), # Separate baseline hazards
                              data = dat_expanded)

  return(cox_long)
}


make_covar_grid <- function(dat) {
  
  #' Make grid of covariate combinations.
  #' 
  #' @inheritParams get_predictor_mats
  #' 
  #' @return Grid of covariate combos.
  #' 
  #' @export
  
  sd_units <- c(0, 1, 2, -1, -2)
  z <- 0 + sd_units * 1 # standard normal
  names(z) <- c("mean", "+1SD", "+2SD","-1SD","-2SD")
  
  if (is.factor(dat$X)) {
    
    x <- c(0, 1)
    names(x) <- c("0", "1")
    
    grid_obj <- expand.grid("X" = (names(x)),
                            "Z" = (names(z)), 
                            stringsAsFactors = F) %>% 
      as.data.frame() %>% 
      mutate(val_X = x[match(X, names(x))],
             val_Z = z[match(Z, names(z))])
                   
    # Maybe unite X and Z?
         
  } else {
    
    x <- 0 + sd_units * 1 # standard normal
    names(x) <- names(z)
    
    grid_obj <- expand.grid("X" = names(x),
                            "Z" = names(z), 
                            stringsAsFactors = F) %>% 
      as.data.frame() %>% 
      tidyr::unite("scen", c(X, Z), sep = "=") %>% 
      dplyr::filter(!(str_detect(.data$scen, "2SD") & 
                        !str_detect(.data$scen, "mean"))) %>% 
      tidyr::separate(.data$scen, c("X", "Z"), sep = "=") %>%            
      dplyr::mutate(val_X = x[match(X, names(x))],
                    val_Z = z[match(Z, names(z))])
    
    # Maybe unite X and Z?
  }
  return(grid_obj)
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
  
  #' @importFrom stats integrate confint
  
  prod <-  function(t) {
    haz_weib(alph_ev, lam_ev, t) * gen_surv_weib(cumhaz_weib(alph_ev, lam_ev, t),
                                                 cumhaz_weib(alph_comp, lam_comp, t))
  }
  
  ci_func <- Vectorize(function(upp) {
    stats::integrate(prod, lower = 0, upper = upp)$value
  })
  
  return(ci_func(t))
}


get_true_cuminc <- function(ev1_pars, 
                            ev2_pars,
                            combo,
                            times) {
  
  #' @title Compute true cumulative incidence
  #' 
  #' @param ev1_pars List arameters for weibull event 1
  #' @param ev2_pars List arameters for weibull event 2
  #' @param combo Covariate combo, e.g. list("val_X" = 1, "val_Z" = 1)
  #' @param times Vector of timepoints to evaluate cuminc at
  #' 
  #' @return True cumulative incidence at times.
  #' 
  #' @export
  
  lam1 <- ev1_pars$h1_0 * exp((ev1_pars$b1 * combo$val_X + 
                                 ev1_pars$gamm1 * combo$val_Z))
  lam2 <- ev2_pars$h2_0 * exp((ev2_pars$b2 * combo$val_X + 
                                 ev2_pars$gamm2 * combo$val_Z))
  
  # only integrate at times!!
  true_pstate2 <- cuminc_weib(alph_ev = ev1_pars$a1, lam_ev = lam1, 
                              alph_comp = ev2_pars$a2, lam_comp = lam2,
                              t = times)
  
  true_pstate3 <- cuminc_weib(alph_ev = ev2_pars$a2, lam_ev = lam2, 
                              alph_comp = ev1_pars$a1, lam_comp = lam1, 
                              t = times)
  
  dato <- cbind.data.frame("times" = times, 
                           "true_pstate2" = true_pstate2,
                           "true_pstate3" = true_pstate3) %>% 
    dplyr::mutate(true_pstate1 = 1 - (true_pstate2 + true_pstate3))
  
  return(dato)
}

