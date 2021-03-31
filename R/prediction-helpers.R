##*****************************##
## Prediction helper functions ##
##*****************************##


# True cumulative incidence calculation -----------------------------------


#' Compute Weibull hazard
#' 
#' @param alph Shape of Weibull distribution
#' @param lam Rate of Weibull distribution
#' @param t Postive scalar or numeric vector, for time
#' 
#' @return Weibull hazard(s) at t as scalar or vector
#' 
#' @export
haz_weib <- Vectorize(function(alph, lam, t) {
  return(alph * lam * t^(alph - 1))
})

#' Compute Weibull cumulative hazard 
#' 
#' @noRd
cumhaz_weib <- Vectorize(function(alph, lam, t) {
  lam * t^alph
})

#' Compute event-free survival (based on Weibull hazards)
#' 
#' @noRd
gen_surv_weib <- Vectorize(function(cumhaz1, cumhaz2) {
  exp(-(cumhaz1 + cumhaz2))
})


#' Compute true cumulative incidence for one of two competing events 
#' 
#' Based on Weibull hazards from both competing events. Integration is
#' done via \code{stats::integrate}.
#' 
#' @param alph_ev Numeric - shape of event of interest
#' @param lam_ev Numeric - rate of event of interest
#' @param alph_comp Numeric - shape of competing event
#' @param lam_comp Numeric - rate of competing event
#' @param t Scalar of vector of (positive) timepoints
#' 
#' @return Scalar, value of cumulative incidence
#' 
#' @export
#' 
#' @examples 
#' time <- seq(0.1, 10, by = 0.1)
#' 
#' # Use baseline parametrisation from sim study
#' cuminc <- cuminc_weib(
#' alph_ev = 0.58,
#' lam_ev = 0.19,
#' alph_comp = 0.53,
#' lam_comp = 0.21,
#' t = time
#' )
#' 
#' plot(
#' x = time, 
#' y = cuminc, 
#' type = "l", 
#' xlab = "Time", 
#' ylab = "Cumulative incidence"
#' )
#' 
cuminc_weib <- function(alph_ev, lam_ev, alph_comp, lam_comp, t) {
  
  prod <-  function(t) {
    haz <- haz_weib(alph = alph_ev, lam = lam_ev, t = t) 
    gen_surv <- gen_surv_weib(
      cumhaz1 = cumhaz_weib(alph = alph_ev, lam = lam_ev, t = t),
      cumhaz2 = cumhaz_weib(alph = alph_comp, lam = lam_comp, t = t)
    )
    haz * gen_surv
  }
  
  ci_func <- Vectorize(function(upp) {
    if (upp == 0) upp <- .Machine$double.eps
    stats::integrate(prod, lower = 0, upper = upp)$value
  })
  
  return(ci_func(t))
}


#' Compute true cumulative incidence (in simulation study)
#' 
#' Relevant only in context of simulation study when generating based on
#' two variables X and Z (also used in vignette)
#' 
#' @param ev1_pars List parameters for weibull event 1
#' @param ev2_pars List parameters for weibull event 2
#' @param combo Covariate combo, e.g. list("val_X" = 1, "val_Z" = 1)
#' @param times Vector of timepoints to evaluate cuminc at
#' 
#' @return True cumulative incidence at given times.
#' 
#' @export
get_true_cuminc <- function(ev1_pars, 
                            ev2_pars,
                            combo,
                            times) {
  
  # Compute rates
  lam1 <- ev1_pars$h1_0 * exp((ev1_pars$b1 * combo$val_X + ev1_pars$gamm1 * combo$val_Z))
  lam2 <- ev2_pars$h2_0 * exp((ev2_pars$b2 * combo$val_X + ev2_pars$gamm2 * combo$val_Z))
  
  # Integrate at times for events 2 and 3
  true_pstate2 <- cuminc_weib(
    alph_ev = ev1_pars$a1, 
    lam_ev = lam1, 
    alph_comp = ev2_pars$a2, 
    lam_comp = lam2,
    t = times
  )
  
  true_pstate3 <- cuminc_weib(
    alph_ev = ev2_pars$a2, 
    lam_ev = lam2, 
    alph_comp = ev1_pars$a1, 
    lam_comp = lam1, 
    t = times
  )
  
  # Combine and return
  dat <- cbind.data.frame(
    "times" = times, 
    "true_pstate2" = true_pstate2, 
    "true_pstate3" = true_pstate3
  ) %>% 
    dplyr::mutate(true_pstate1 = 1 - (true_pstate2 + true_pstate3))
  
  return(dat)
}


# Predicting for covariate combinations -----------------------------------


#' Make grid of covariate combinations (i.e. reference patients to predict)
#' 
#' @inheritParams get_predictor_mats
#' 
#' @return Grid of covariate combos.
#' 
#' @noRd
make_covar_grid <- function(dat) {
  
  X <- Z <- NULL
  
  sd_units <- c(0, 1, 2, -1, -2)
  z <- 0 + sd_units * 1 # standard normal
  names(z) <- c("mean", "+1SD", "+2SD","-1SD","-2SD")
  
  if (is.factor(dat$X)) {
    
    x <- c(0, 1)
    names(x) <- c("0", "1")
    
    grid_obj <- expand.grid(
      "X" = (names(x)),
      "Z" = (names(z)), 
      stringsAsFactors = FALSE
    ) %>% 
      as.data.frame() %>% 
      dplyr::filter(!(stringr::str_detect(.data$Z, "2SD"))) %>% 
      dplyr::mutate(
        val_X = x[match(X, names(x))],
        val_Z = z[match(Z, names(z))]
      )
    
  } else {
    
    x <- 0 + sd_units * 1 # standard normal
    names(x) <- names(z)
    
    grid_obj <- expand.grid(
      "X" = names(x),
      "Z" = names(z), 
      stringsAsFactors = FALSE
    ) %>% 
      as.data.frame() %>% 
      tidyr::unite("scen", c(X, Z), sep = "=") %>% 
      dplyr::filter(
        !(stringr::str_detect(.data$scen, "2SD")),
        !(stringr::str_detect(.data$scen, "mean"))
      ) %>% 
      dplyr::add_row(scen = "mean=mean") %>% 
      tidyr::separate(.data$scen, c("X", "Z"), sep = "=") %>%            
      dplyr::mutate(
        val_X = x[match(X, names(x))],
        val_Z = z[match(Z, names(z))]
      )
    
  }
  return(grid_obj)
}

#' Compute state probabilities based on mstate competing risks model
#' 
#' This is a dplyr-based version of `mstate::probtrans`, without computing
#' standard errors. Yields exactly the same results.
#' 
#' @noRd
get_state_probs <- function(obj, 
                            combo,
                            cox_long) {
  
  # For checks
  trans <- Haz <- `1` <- `2` <- H_REL <- H_NRM <- haz_REL <- haz_NRM <- NULL
  hazsum <- time <- pstate1 <- pstate2 <- pstate3 <- EFS_min1 <- NULL
  
  # Read in baseline cumulative hazard
  df_baseHaz <- obj$Haz
  
  # Compute linear predictors
  lp_REL <- as.numeric(cox_long$coefficients[1:2] %*% c(combo$val_X, combo$val_Z))
  lp_NRM <- as.numeric(cox_long$coefficients[3:4] %*% c(combo$val_X, combo$val_Z))
  
  # Start computation of probabilities
  df_cumincs <- df_baseHaz %>% 
    tidyr::pivot_wider(names_from = trans, values_from = Haz) %>% 
    dplyr::rename("H_REL" = `1`, "H_NRM" = `2`) %>% 
    
    # Multiply estimated cumulative hazards by Linear predictor
    dplyr::mutate(H_REL = H_REL * exp(lp_REL), H_NRM = H_NRM * exp(lp_NRM)) %>% 
    
    # Compute hazards and EFS
    dplyr::mutate(
      haz_REL = diff(c(0, H_REL)),
      haz_NRM = diff(c(0, H_NRM)),
      hazsum = haz_REL + haz_NRM,
      pstate1 = cumprod(1 - hazsum)
    ) %>% 
    
    # Add time-point zero
    dplyr::add_row(time = 0, haz_REL = 0, haz_NRM = 0, pstate1 = 1) %>% 
    dplyr::arrange(time) %>% 
    
    # Compute cumulative incidences 
    dplyr::mutate(
      EFS_min1 = c(1, pstate1[-length(pstate1)]),
      pstate2 = cumsum(haz_REL * EFS_min1),  
      pstate3 = cumsum(haz_NRM * EFS_min1)
    ) %>% 
    
    # Keep essentials
    dplyr::select(time, pstate1, pstate2, pstate3)
  
  return(df_cumincs)
}


#' Predicting cumulative incidence for grid of patients
#' 
#' (Simulation study helper)
#' 
#' @inheritParams get_true_cuminc
#' @param cox_long Cox model returned by setup_mstate
#' @param times Prediction horizon(s) eg. c(2, 5, 10)
#' @param grid_obj Grid object
#' 
#' @return Df of predictions
#' 
#' @noRd
get_preds_grid <- function(cox_long,
                           grid_obj,
                           times,
                           ev1_pars,
                           ev2_pars) {
  
  # Set-up for baseline (covariate values = 0)
  tmat <- mstate::trans.comprisk(K = 2, names = c("Rel", "NRM"))

  baseline_dat <- data.frame(
    "X.1" = c(0, 0),
    "X.2" = c(0, 0), 
    "Z.1" = c(0, 0), 
    "Z.2" = c(0, 0), 
    "trans" = c(1, 2), 
    "strata" = c(1, 2)
  )
  
  # Run msfit once to get baseline cumulative hazards for all causes                        
  msfit_baseline <- mstate::msfit(
    object = cox_long,
    newdata = baseline_dat,
    trans = tmat
  )
  
  # Get predicted state probabilites for all covariate combos
  preds_grid <- purrr::map_dfr(1:nrow(grid_obj), function(row) {
    
    combo <- grid_obj[row, ]
    
    # Get state probabilities at all time points, for this combo
    preds_test <- get_state_probs(
      msfit_baseline, 
      combo = combo, 
      cox_long = cox_long
    )
    
    # Summarise for prediction horizons of interest
    summ <- purrr::map_dfr(
      times, 
      ~ {
        preds_test %>% 
          dplyr::filter(time <= .x) %>% 
          dplyr::slice(dplyr::n()) %>% 
          dplyr::select(times = time, pstate1, pstate2, pstate3) %>% 
          dplyr::mutate(times = .x) %>% 
          as.data.frame()
      } 
    )
    
    # Compute true cumulative incidence at those horizons
    true_CI <- get_true_cuminc(
      ev1_pars = ev1_pars, 
      ev2_pars = ev2_pars,
      combo = combo, 
      times = times
    )
    
    # Join the true and predicted state probabilities
    res <- cbind.data.frame(summ, "X" = combo$X, "Z" = combo$Z) %>% 
      dplyr::left_join(true_CI, by = "times") %>% 
      dplyr::mutate_if(is.factor, as.character) # avoid bind_rows warnings
    
    return(res)
  })
  
  return(preds_grid)
}
