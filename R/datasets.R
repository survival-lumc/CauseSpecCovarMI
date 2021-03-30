#' Summarised simulation results - regression coefficients
#'
#' @docType data
#' 
#' @format 
#' \describe{
#' \item{var}{Variable for which a regression coefficient is being estimated, for example
#' "X.1" corresponds to the cause-specific effect of X on event 1 (Relapse)}
#' \item{m}{Number of imputed datasets}
#' \item{analy}{Method indicator, please refer to manuscript for abbreviations used}
#' \item{true}{True data-generating value for regression coefficient}
#' \item{scen_num}{Number label for scenario}
#' \item{n}{Sample size in simulated datasets}
#' \item{est}{Estimated regression coefficient}
#' \item{se}{Estimate model standard error (SE) for regression coefficient}
#' \item{se_mcse}{Monte carlo standard error for model SE}
#' \item{emp_se}{Empirical SE}
#' \item{cover}{Coverage}
#' \item{bias}{Bias}
#' \item{rmse}{Coverage}
#' \item{rmse_mcse}{Coverage}
#' \item{warns}{Number of warnings across all simulations of a particular scenario, for 
#' a coefficient. Generally corresponds to number of smcfcs rejection sampling failures.}
#' \item{bias_mcse}{Monte carlo error of bias}
#' \item{cover_mcse}{Monte carlo error of coverage}
#' \item{prop_miss}{Proportion of missing values in X}
#' \item{beta1}{Data-generating ffect of X on event 1}
#' \item{miss_mech}{Missingness mechanism}
#' \item{X_level}{Measurement level of X}
#' \item{rho}{Correlation between X and Z}
#' \item{eta1}{Strength of missingness mechanism for non-MCAR scenarios}
#' \item{haz_shape}{Shapes of baseline hazard for competing events, either "similar"
#'  or "different"}
#' }
#'
#' @usage regr_results
#'
#' @keywords datasets
"regr_results"

#' Summarised simulation results - predictions
#'
#' @docType data
#' 
#' @format 
#' \describe{
#' \item{analy}{Method indicator, please refer to manuscript for abbreviations used}
#' \item{m}{Number of imputed datasets}
#' \item{combo-X_Z}{Covariate combination defining reference patient for which
#' we predict}
#' \item{times}{Prediction horizon (in years)}
#' \item{state}{State for which we estimate probability at times}
#' \item{true}{True data-generating cumulative incidence at a horizon, and for 
#' a certain reference patient}
#' \item{scen_num}{Number label for scenario}
#' \item{n}{Sample size in simulated datasets}
#' \item{prob}{Estimated probability}
#' \item{emp_se}{Empirical SE of the probabilities}
#' \item{bias}{Bias}
#' \item{rmse}{Coverage}
#' \item{rmse_mcse}{Coverage}
#' \item{bias_mcse}{Monte carlo error of bias}
#' \item{prop_miss}{Proportion of missing values in X}
#' \item{beta1}{Data-generating ffect of X on event 1}
#' \item{miss_mech}{Missingness mechanism}
#' \item{X_level}{Measurement level of X}
#' \item{rho}{Correlation between X and Z}
#' \item{eta1}{Strength of missingness mechanism for non-MCAR scenarios}
#' \item{haz_shape}{Shapes of baseline hazard for competing events, either "similar"
#'  or "different"}
#' }
#'
#' @usage preds_results
#'
#' @keywords datasets
"preds_results"
