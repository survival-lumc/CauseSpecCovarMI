##***************************##
## General utility functions ##
##***************************##


# Imports from other packages ---------------------------------------------


#' @importFrom magrittr `%$%`
#' @importFrom rlang .data
#' @importFrom data.table .N .I ':=' .SD
NULL


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @export
NULL


#' Imported \code{survival::Surv}
#'
#' See \code{survival::\link[survival:Surv]{\%>\%}} for details.
#'
#' @name Surv
#' @rdname Surv
#' @keywords internal
#' @importFrom survival Surv
#' @export
NULL


#' Imported \code{survival::strata}
#'
#' See \code{survival::\link[survival:strata]{\%>\%}} for details.
#'
#' @name strata
#' @rdname strata
#' @keywords internal
#' @importFrom survival strata
#' @export
NULL


#' Imported pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @export
NULL


# General utilities -------------------------------------------------------


#' Silence function printing
#' 
#' Takes a function/expression that by default uses print() or cat(), and stops it. 
#' We use it for \code{MIcombine()} or \code{smcfcs()}.
#' 
#' @source \link{http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html}
#' 
#' @param expr Expression to silence
#' 
#' @noRd
quiet <- function(expr) { 
  
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(expr)) 
} 


#' Record warnings
#' 
#' Used to capture number of rejection sampling failures from \code{smcfcs}.
#' 
#' @param expr Expression from which to record possible warning
#' 
#' @return List with expression and associated warning
#' 
#' @noRd
record_warning <- function(expr) {
  
  warn <- NULL
  
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
#' Used to collapse scenario identifiers in one column.
#' 
#' @export
#' @noRd
add_scen_details <- function(scenario,
                             seed,
                             rep_num) {
  
  name <- NULL
  
  scen_dat <- data.frame(t(scenario)) %>% 
    tibble::rownames_to_column(var = "name") %>% 
    dplyr::filter(!(name %in% c("pilot", "seed"))) %>% 
    tidyr::unite("scen", c(1, 2), sep = "=") 
  
  scen_collapse <- paste(scen_dat$scen, collapse = "-")
  rep <- paste0("rep=", rep_num)
  seed <- paste0("seed=", seed)
  
  return(paste(c(scen_collapse, rep, seed), collapse = "-"))
}


#' Format CCA and reference analyses
#' 
#' Essentially a version \code{broom::tidy} - made before I knew of the existence
#' of the broom package.
#' 
#' @param mod Cause-specific Cox model as returned by çode{survival::coxph}
#' @inheritParams pool_diffm
#' 
#' @return Formatted dataframe with model summary
#'
#' @noRd
summarise_ref_CCA <- function(mod,
                              analy) {
  
  summ_ref_CCA <- cbind(summary(mod)$coefficients, stats::confint(mod)) %>% 
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
#' Similar to \code{broom::tidy} - but specifically for the probabilities
#' returned by \code{mstate::probtrans}. Only works for competing risks data.
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
    dplyr::mutate(
      state = dplyr::case_when(
        stringr::str_detect(.data$state, "pstate1_true_pstate1") ~ "1",
        stringr::str_detect(.data$state, "pstate2_true_pstate2") ~ "2",
        stringr::str_detect(.data$state, "pstate3_true_pstate3") ~ "3"
      )
    ) %>% 
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


# Adapted from calc_absolute - without dplyr::pull for speed,
# operating solely on vectors
# 
# 

#' Approximate jackknife estimate for MCSE of RMSE
#' 
#' Adapted from the simhelpers CRAN package, see \link{https://cran.r-project.org/web/packages/simhelpers/vignettes/MCSE.html}.
#' This version takes vectors as inputs instead.
#' 
#' @source \link{https://github.com/meghapsimatrix/simhelpers/blob/master/R/calc_absolute.R}
#' 
#' @export
rmse_mcse <- function(estimates, true, K) {
  
  # Keep first true value
  true_param <- true[1]
  
  # Calculate elements of rmse mcse
  t_bar <- mean(estimates) 
  var_t <- stats::var(estimates) 
  t_bar_j <- (1 / (K - 1)) * (K * t_bar - estimates) 
  bias_j_sq <- (t_bar_j - true_param)^2 
  s_sq_t_j <- (1 / (K - 2)) * ((K - 1) * var_t - (K / (K - 1)) * (estimates - t_bar)^2) 
  rmse_j <- sqrt(bias_j_sq + s_sq_t_j) 
  mse <- mean((estimates - true_param)^2) 
  
  # Calculate rmse and mcse
  rmse <- sqrt(mse)
  mcse <- sqrt(((K - 1) / (K)) * sum((rmse_j - rmse)^2))
  
  return(mcse)
}


#' Helper function to extract scenario number
#' 
#' @export
#' @noRd
extract_scen_num <- function(dat) {
  
  scen_summary <- NULL
  
  dat[, ':=' (
    scen_num = gsub(
      pattern = ".*(scen_num=)|(-rep).*$", 
      replacement = "", 
      x = scen_summary
    ),
    scen_summary = gsub(
      pattern = "(-scen_num).*$", 
      replacement = "", 
      x = scen_summary
    )
  )]
  
  return(dat)
}

#' Helper function for formatting full simulation results before summarising
#' 
#' @export
#' @noRd
format_scen_summary <- function(dat) {
  
  # For checks
  analy <- prop_miss <- haz_shape <- beta1 <- eta1 <- NULL
  scen_summary <- miss_mech <- X_level <- rho <- n <- NULL
  
  # Labels for variables in the "scen_summary" column
  labs_scens <- c(
    "n", 
    "prop_miss", 
    "beta1", 
    "miss_mech", 
    "X_level",
    "rho", 
    "eta1", 
    "haz_shape"
  )
  
  # Adjust for eta minus label
  dat[, scen_summary := gsub(
    pattern = "eta1=-", 
    replacement = "eta1=min", 
    x = scen_summary
  )]
  
  # Separate the scen_summary column
  dat[, (labs_scens) := data.table::tstrsplit(scen_summary, split = "-")] 
  
  # Relevel new variables - ready for analysis
  dat[, ':=' (
    analy = factor(
      analy, levels = c("ref", "CCA", "ch1", "ch12", "ch12_int", "smcfcs")
    ),
    prop_miss = factor(
      prop_miss, 
      levels = c("prop_miss=0.1","prop_miss=0.5"),
      labels = c("10%", "50%")
    ),
    haz_shape = factor(
      haz_shape,
      levels = c("haz_shape=similar", "haz_shape=different"),
      labels = c("similar", "different")
    ),
    beta1 = factor(
      beta1,
      levels = c("beta1=0", "beta1=0.5", "beta1=1"),
      labels = c("0", "0.5", "1")
    ),
    eta1 = factor(
      eta1,
      levels = c("eta1=NA", "eta1=min1", "eta1=min2"),
      labels = c("None", "Weak", "Strong")
    ),
    miss_mech = factor(
      miss_mech,
      levels = c("miss_mech=MCAR", "miss_mech=MAR", 
                 "miss_mech=MAR_GEN", "miss_mech=MNAR"),
      labels = c("MCAR", "MAR", "MAR_GEN", "MNAR") # change to MAR-T ?
    ),
    X_level = factor(
      X_level,
      levels = c("X_level=continous", "X_level=binary"),
      labels = c("continuous", "binary")
    ),
    rho = factor(
      rho,
      levels = "rho=0.5",
      labels = "0.5"
    ),
    n = factor(
      n,
      levels = c("n=500", "n=2000"),
      labels = c("500", "2000")
    ),
    
    # Remove original summary variable
    scen_summary = NULL
  )]
  
  return(dat)
}
