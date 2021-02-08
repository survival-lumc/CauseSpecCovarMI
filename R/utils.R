##***************************##
## General utility functions ##
##***************************##


# General -----------------------------------------------------------------


#' Silence function printing
#' 
#' Takes a function/expression that by default uses 
#' print() or cat(), and stops it. We use ot for MIcombine() or smcfcs()
#' 
#' Source: http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
#' 
#' @param expr Expression to silence
#' 
#' @noRd
quiet <- function(expr) { 
  
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(expr)) 
} 


#' Record rejection sampling failure smcfcs
#' 
#' @param expr Expression from which to record possible error
#' 
#' @return List with expression and associated warning
#' 
#' @noRd
record_warning <- function(expr) {
  
  value <- withCallingHandlers(expr, warning = function(w) {
    warn <<- w
    invokeRestart("muffleWarning")
  })
  
  # Check if there was a warning at all
  warn_val <- ifelse(exists("warn"), as.character(warn), as.character(0))
  
  return(list(value = value, warning = warn_val))
}


# Formatting simulation results -------------------------------------------


#' Appending scenario details as column
#' 
#' @noRd
add_scen_details <- function(scenario,
                             seed,
                             rep_num) {
  
  scen_dat <- data.frame(t(scenario)) %>% 
    tibble::rownames_to_column(var = "name") %>% 
    dplyr::filter(!(name %in% c("pilot", "seed"))) %>% 
    tidyr::unite("scen", 1:2, sep = "=") 
  
  scen_collapse <- paste(scen_dat$scen, collapse = "-")
  rep <- paste0("rep=", rep_num)
  seed <- paste0("seed=", seed)
  
  return(paste(c(scen_collapse, rep, seed), collapse = "-"))
}


#' Format CCA and ref analyses
#' 
#' @param mod Cause-specific model
#' @inheritParams pool_diffm
#' 
#' @return Formatted df
#' 
#' @noRd
summarise_ref_CCA <- function(mod,
                              analy) {
  
  summ_ref_CCA <- cbind(summary(mod)$coefficients,
                        stats::confint(mod)) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("var") %>% 
    dplyr::rename(
      estimate = .data$coef, 
      std.error = .data$`se(coef)`,
      p.value = .data$`Pr(>|z|)`
    ) %>% 
    dplyr::select(-.data$`exp(coef)`, -.data$z) %>% 
    dplyr::mutate(m = 0, analy = analy)
  
  return(summ_ref_CCA)
}


#' Format predictions of CCA and ref
#' 
#' @noRd
preds_CCA_ref <- function(preds,
                          analy) {
  
  res <- preds %>% 
    tidyr::pivot_longer(
      cols = .data$pstate1:.data$pstate3, 
      names_to = "state_est", 
      values_to = "prob"
    ) %>% 
    tidyr::pivot_longer(
      cols = .data$true_pstate2:.data$true_pstate1, 
      names_to = "state_true", 
      values_to = "true"
    ) %>% 
    
    # Make single state variable
    tidyr::unite("state", .data$state_est, .data$state_true) %>% 
    dplyr::mutate(state = dplyr::case_when(
      stringr::str_detect(.data$state, "pstate1_true_pstate1") ~ "1",
      stringr::str_detect(.data$state, "pstate2_true_pstate2") ~ "2",
      stringr::str_detect(.data$state, "pstate3_true_pstate3") ~ "3"
    )) %>% 
    dplyr::filter(!is.na(.data$state)) %>% 
    tidyr::unite("combo-X_Z", .data$X, .data$Z, sep = "_X-Z_") %>% 
    
    # Labels
    dplyr::rename(p_pool = "prob") %>% 
    dplyr::mutate(
      analy = analy,
      m = 0,
      sq_err = (.data$p_pool - .data$true)^2
    )
  
  return(res)
}
