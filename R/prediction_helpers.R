##*****************************##
## Prediction helper functions ##
##*****************************##


get_state_probs <- function(obj, # needs to be baseline (covars 0)
                            combo,
                            cox_long) {
  
  # Read in baseline cumulative hazard
  df_baseHaz <- obj$Haz
  
  # Compute linear predictors
  lp_REL <- as.numeric(
    cox_long$coefficients[1:2] %*% c(combo$val_X, combo$val_Z)
  )
  
  lp_NRM <- as.numeric(
    cox_long$coefficients[3:4] %*% c(combo$val_X, combo$val_Z)
  )
  
  # CumInc pipeline:
  df_cumincs <- df_baseHaz %>% 
    tidyr::pivot_wider(names_from = trans, 
                       values_from = Haz) %>% 
    dplyr::rename("H_REL" = `1`,
                  "H_NRM" = `2`) %>% 
    
    # Multiply estimated cumulative hazards by Linear predictor
    dplyr::mutate(H_REL = H_REL * exp(lp_REL),
                  H_NRM = H_NRM * exp(lp_NRM)) %>% 
    
    # Compute hazards and EFS
    dplyr::mutate(haz_REL = diff(c(0, H_REL)),
                  haz_NRM = diff(c(0, H_NRM)),
                  hazsum = haz_REL + haz_NRM,
                  pstate1 = cumprod(1 - hazsum)) %>% 
    
    # Add time-point zero
    dplyr::add_row(time = 0, haz_REL = 0, haz_NRM = 0, pstate1 = 1) %>% 
    dplyr::arrange(time) %>% 
    
    # Compute cumulative incidences 
    dplyr::mutate(EFS_min1 = c(1, pstate1[-length(pstate1)]),
                  pstate2 = cumsum(haz_REL * EFS_min1),  
                  pstate3 = cumsum(haz_NRM * EFS_min1)) %>% 
    
    # Keep essentials
    dplyr::select(time, pstate1, pstate2, pstate3)
  
  return(df_cumincs)
}



get_preds_grid <- function(cox_long,
                           grid_obj,
                           times,
                           ev1_pars,
                           ev2_pars) {
  
  #' Predicting grids
  #' 
  #' @param cox_long Cox model returned by setup_mstate
  #' @param times Prediction horizon(s) eg. c(2, 5, 10)
  #' @inheritParams get_true_cuminc
  #' @param grid_obj Grid object
  #' 
  #' @return Df of predictions
  #' 
  #' @export
  
  # Set-up for baseline (covariate values = 0)
  tmat <- mstate::trans.comprisk(2, c("Rel", "NRM"))

  baseline_dat <- data.frame(X.1 = c(0, 0),
                             X.2 = c(0, 0), 
                             Z.1 = c(0, 0), 
                             Z.2 = c(0, 0), 
                             trans = c(1, 2), 
                             strata = c(1, 2))
  
  # Run msfit once to get baseline cumulative hazards for all causes                        
  msfit_baseline <- mstate::msfit(cox_long,
                                  newdata = baseline_dat,
                                  trans = tmat)
                         
  
  # Get predicted state probabilites for all covariate combos
  preds_grid <- purrr::map_dfr(1:nrow(grid_obj), function(row) {
    
    combo <- grid_obj[row, ]
    
    # Get state probabilities at all time points, for this combo
    preds_test <-  SimsCauseSpecCovarMiss::get_state_probs(
      msfit_baseline, 
      combo = combo, 
      cox_long = cox_long
    )
    
    # Summarise for prediction horizions of interest
    summ <- purrr::map_dfr(
      times, 
      ~ preds_test %>% 
        dplyr::filter(time <= .x) %>% 
        dplyr::slice(dplyr::n()) %>% 
        dplyr::select(times = time, pstate1, pstate2, pstate3) %>% 
        dplyr::mutate(times = .x) %>% 
        as.data.frame()
    )
    
    # Compute true cumulative incidence at those horizons
    true_CI <- SimsCauseSpecCovarMiss::get_true_cuminc(
      ev1_pars = ev1_pars, 
      ev2_pars = ev2_pars,
      combo = combo, 
      times = times
    )
    
    # Join the true and predicted state probabilities
    res <- cbind.data.frame(summ,
                            "X" = combo$X, 
                            "Z" = combo$Z) %>% 
      dplyr::left_join(true_CI, by = "times") %>% 
      dplyr::mutate_if(is.factor, as.character) # avoid bind_rows warnings
    
    return(res)
  })
  
  return(preds_grid)
}
