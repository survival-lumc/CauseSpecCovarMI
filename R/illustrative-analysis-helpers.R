##****************************************##
## Helper functions illustrative analysis ##
##****************************************##



#' Extract variable names from right side of model formula
#' 
#' @export
#' @noRd
extract_rhs_varnames <- function(form, dat) {
  
  # Extract rhs terms
  coef_names <- attr(x = terms(form), which = "term.labels")
  
  # Get column names to match
  colnames_pattern <- paste0("(", paste(colnames(dat), collapse = "|"), ")")
  
  # Match
  rhs_varnames <- unique(
    regmatches(x = coef_names, m = regexpr(text = coef_names, pattern = colnames_pattern))
  )
  
  return(rhs_varnames)
}

#' Choose standard reference patient to predict
#' 
#' Based on either most common or reference factor levels in sample for
#' categorical variable, or median/mean for continuous variables.
#' 
#' @export
choose_standard_refpat <- function(col,
                                   contin_action = c("median", "mean"),
                                   categ_action = c("most_common", "reference")) {
  
  # Convert to factor if characters
  if (is.character(col)) col <- factor(col)
  
  # Match args
  contin <- match.arg(contin_action)
  categ <- match.arg(contin_action)
  
  # Set if continuous
  if (is.numeric(col)) {
    mean_val <- mean(col, na.rm = TRUE)
    median_val <- median(col, na.rm = TRUE)
    val <- ifelse(contin == "mean", mean_val, median_val) 
  } else {
    common_cat <- names(which.max(table(col, useNA = "no")))
    reference_cat <- levels(col)[1]
    val <- ifelse(categ == "most_common", common_cat, reference_cat)
    val <- factor(val, levels = levels(col))
  }
  
  return(val)
}

#' Prepare reference patient to predict with probtrans
#' 
#' Utility function which allows to prepare one row of data.frame
#' to be fed into \code{mstate::msfit} for prediction.
#' 
#' @param refpat A single-row data.frane containing the covariate
#' values to predict
#' @param tmat Transition matrix
#' @param covs Covariates used in the cause-specific Cox models
#' 
#' @export
make_mstate_refpat <- function(refpat, tmat, covs) {
  
  # Get number of transitioins
  n_trans <- max(tmat, na.rm = TRUE)
  
  # Copy ref pat n_trans times (one per transition)
  refpat_new <- do.call("rbind", replicate(n_trans, refpat, simplify = FALSE))
  
  # Add trans variable and set attributes
  refpat_new$trans <- 1:n_trans
  attr(refpat_new, "trans") <- tmat
  class(refpat_new) <- c("msdata", "data.frame")
  
  # Expand now
  refpat_expanded <- mstate::expand.covs(refpat_new, covs, longnames = FALSE)
  refpat_expanded$strata <- 1:n_trans
  return(refpat_expanded)
}

#' Utility to add reference factor levels
#' 
#' @export
#' @noRd
reflevels_add_summary <- function(summ, dat, form, term_col = "term") {
  
  # Get predictors
  preds <- attr(terms(form), "term.labels")
  
  # Identify factors
  ref_levels <- dat[, lapply(.SD, function(col) levels(col)[1]), .SDcols = is.factor] %>% 
    data.table::transpose(keep.names = "variable") %>% 
    data.table::setnames(old = "V1", new = "coef")
  
  ref_levels[, term := paste0(variable, coef)]
  ref_levels[, setdiff(names(ref_levels), "term") := NULL]
  
  return(rbind(data.table::data.table(summ), ref_levels, fill = TRUE))
}


# Pooling predictions -----------------------------------------------------

#' Helper function to run models in the illustrative analysis
#' 
#' @export
#' @noRd
run_mds_model <- function(form,
                          tmat,
                          dat) {
                              
  
  # Get predictor names from formula
  predictors <- extract_rhs_varnames(form, dat)
  
  # Prepare impdat for mstate model
  if (!(any(class(dat) %in% "data.table"))) dat <- data.table::data.table(dat)
  dat[, ':=' (ev1 = ci_s_allo1 == 1, ev2 = ci_s_allo1 == 2)]
  dat <- as.data.frame(dat)
  
  dat_msprepped <- mstate::msprep(
    time = c(NA, "ci_allo1", "ci_allo1"),
    status = c(NA, "ev1", "ev2"), 
    data = dat,
    trans = tmat,
    keep = predictors
  ) 
  
  dat_expanded <- mstate::expand.covs(
    dat_msprepped, predictors, append = TRUE, longnames = FALSE
  )
  
  # Run do.call due to non-standard formula eval
  #mod <- do.call(survival::coxph, list("formula" = form, "data" = dat_expanded))
  mod <- survival::coxph(formula = form, data = dat_expanded, model = TRUE)
  
  return(mod)
}


#' Helper function to obtain predictions in the illustrative analysis
#' 
#' (For each imputed dataset)
#' 
#' @export
#' @noRd
predict_mds_model <- function(mod,
                              ref_pats,
                              tmat, 
                              horizon) {
  
  # Iterate over refpaths
  probs_df <- purrr::map_dfr(
    .x = ref_pats,
    .f = ~ {
      
      # First with one refpat
      msf_obj <- mstate::msfit(mod, newdata = .x, trans = tmat)
      
      # Run probtrans
      pt_obj <- mstate::probtrans(msf_obj, predt = 0)[[1]]
      
      # Get probability at horizon
      probs <- data.table::last(pt_obj[pt_obj$time <= horizon, ])
      
      # res
      res <- cbind.data.frame(
        "state" = 1:3,
        "prob" = with(probs, c(pstate1, pstate2, pstate3)),
        "se" = with(probs, c(se1, se2, se3)),
        "horiz" = horizon
      )
      
      return(res)
    }, .id = "ref_pat"
  )
  
  return(probs_df)
}


#' @export
#' @noRd
cloglog <- Vectorize(function(x) log(-log(1 - x)))

#' @export
#' @noRd
inv_cloglog <- Vectorize(function(x) 1 - exp(-exp(x)))


#' Pooling probabilities based on complementary log-log transformation
#' 
#' Method described in paper of morisot and colleagues
#' 
#' @param preds_list A list of length equal to number of imputed datasets,
#' containing the imputation-specific predictions. Each element should be a dataframe
#' containing columns "prob" (probability), "se" (standard error of probability) and 
#' any other variables which identify groups of predictions (to be used in by_vars)
#' @param by_vars Vector of variable names to pool across
#' 
#' @export
pool_morisot <- function(preds_list, # add p_var and se_var, also confint
                         by_vars) {
  
  # Bind the predictions  
  preds_full <- data.table::rbindlist(preds_list)

  # Step 1
  preds_full[, ':=' (
    p  = prob,
    p_trans = cloglog(prob),
    var_p = (se)^2
  )] 
  
  preds_full[, Ui := var_p / (log(1 - p) * (1 - p))^2]
  
  preds_summ <- preds_full[, .(
    m = .N,
    Ubar = mean(Ui),
    B = var(p_trans),
    Qbar = mean(p_trans)
  ), by = by_vars]
  
  preds_summ[, total_var := Ubar + (1 + m^-1) * B, by = by_vars]
  
  preds_final <- preds_summ[, .(
    p_pooled = inv_cloglog(Qbar),
    CI_low = inv_cloglog(Qbar - qnorm(0.975) * sqrt(total_var)),
    CI_upp = inv_cloglog(Qbar + qnorm(0.975) * sqrt(total_var))
  ), by = by_vars]
  
  return(preds_final)
}