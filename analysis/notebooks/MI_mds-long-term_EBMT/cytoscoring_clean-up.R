##***********************************##
##  Cytoscoring cleaned-up analysis  ##
##***********************************##

# Load libraries
pacman::p_load(
  sjlabelled,
  janitor,
  data.table,
  magrittr,
  purrr,
  survival, 
  mstate,
  mice, 
  smcfcs,
  mitools,
  naniar,
  ggpubr, # will load ggplot with it
  gtsummary
)


# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2)) 

# Reprod.
set.seed(1984)


# Read-in data ------------------------------------------------------------



# Read-in data with extra vars
dat_orig <- sjlabelled::read_spss(
  path = "analysis/data/raw_data/sAML_final_20160630_withcytoscoring_b.sav" 
) %>% 
  
  # Clean variable names (i.e. bring lower case, with _ in between)
  janitor::clean_names() %>% 
  
  # Make into data.table
  data.table()


# Define variables to keep - from paper
vars_keep <- c(
  
  # Outcome vars
  "ci_s", 
  "ci_t", 
  "srv_s", 
  "srv_t",
  
  # aa_auto and covariates
  "aa_auto", 
  "age_allo1", 
  "classification_diag1",
  "donor", 
  "source", 
  "tcelexvivo_allo1", 
  "cytog_recat",
  "stage_cr_allo1",
  "karnofsk_allo1", 
  "cmvpat",
  "ric"
)


# Prepare data ------------------------------------------------------------


# Factors that will use their SPSS value labels 
vars_val_labs <- c(
  "classification_diag1",
  "donor", 
  "source", 
  "tcelexvivo_allo1", 
  "cytog_recat",
  "stage_cr_allo1",
  "cmvpat",
  "ric"
)

# Prepare
dat_reg <- data.table::copy(dat_orig) %>% 
  
  # Select variables
  .[, ..vars_keep] %>% 
  
  # Set time in years
  .[, ':=' (
    ci_t = ci_t / 12,
    srv_t = srv_t / 12
  )] %>% 
  
  # Artifical censor at 3y
  .[ci_t >= 3, ':=' (
    ci_t = 3,
    ci_s = 0
  )] %>% 
  
  # Artifical censor at 3y
  .[srv_t >= 3, ':=' (
    srv_t = 3,
    srv_s = 0
  )] %>% 
  
  # For some vars, use the spss values labels as factor labels
  .[, (vars_val_labs) := lapply(.SD, sjlabelled::as_label), 
    .SDcols = vars_val_labs] %>% 
  
  # Cytogenetics missing code 
  .[cytog_recat == "Missing", cytog_recat := NA_character_] %>% 
  
  # Make cytogenetics orders
  .[, cytog_recat := as.ordered(cytog_recat)] %>% 
  
  # Dichtomise karnof at 90
  .[, karnofsk_allo1 := cut(
    x = karnofsk_allo1, 
    breaks = c(0, 80, 100),
    include.lowest = TRUE, 
    labels = c("<=80", ">80") 
  )] %>% 
  
  # Reorder Karnofsky
  .[, karnofsk_allo1 := relevel(
    karnofsk_allo1, ref = ">80"
  )] %>% 
  
  # Age in decades
  .[, age_allo1 := age_allo1 / 10]  %>% 

  # Remove 135 patients with missing ci_s info
  .[!is.na(ci_s),]


# Drop unused factor levels 
which_factors <- names(dat_reg)[sapply(dat_reg, is.factor)]
dat_reg[, (which_factors) := lapply(.SD, droplevels), .SDcols = which_factors]  


# Conserve var labels from original data (for mutated vars)
dat_reg <- data.table::copy(
  sjlabelled::copy_labels(df_new = dat_reg, df_origin = dat_orig)
)  



# CCA ---------------------------------------------------------------------


predictors <- c("classification_diag1",
                "age_allo1", 
                "stage_cr_allo1",
                "donor", 
                "source", 
                "cmvpat",
                "ric",
                "karnofsk_allo1", 
                "cytog_recat",
                "tcelexvivo_allo1")

names(predictors) <- predictors

# RHS of formula
preds <- paste(predictors, collapse = " + ")
preds 

# For relapse
form_rel <- as.formula(paste0("Surv(ci_t, ci_s == 1) ~ ", preds))
form_rel


mod_rel <- coxph(formula = form_rel, data = dat_reg)

broom::tidy(mod_rel) %>%
  broom::fix_data_frame() %>% 
  dplyr::select(term, estimate, conf.low, conf.high, p.value) %>% 
  dplyr::mutate(estimate = exp(estimate))


# For nrm
form_nrm <- as.formula(paste0("Surv(ci_t, ci_s == 2) ~ ", preds))

mod_nrm <- coxph(formula = form_nrm, data = dat_reg)
mod_nrm <- update(mod_nrm, . ~ . - ric)

broom::tidy(mod_nrm) %>%
  broom::fix_data_frame() %>% 
  dplyr::select(term, estimate, conf.low, conf.high, p.value) %>% 
  dplyr::mutate(estimate = exp(estimate))



# Prep MI -----------------------------------------------------------------


# First, for competing risks, we make indicators for comp events and
# compute cumulative hazards

# Get event indictors
dat_reg[, ':=' (
  ev1 = ifelse(ci_s == 1, 1, 0),
  ev2 = ifelse(ci_s == 2, 1, 0)
)]

# Compute cumulative hazards
dat_reg[, ':=' (
  H1 = mice::nelsonaalen(data = dat_reg, timevar = ci_t, statusvar = ev1),
  H2 = mice::nelsonaalen(data = dat_reg, timevar = ci_t, statusvar = ev2)
)]


# Which vars actually have some data missing?
var_names_miss <- naniar::miss_var_which(dat_reg)

# Set methods accordingly
meth_miss <- setNames(
  object = c("logreg",  "polr", rep("logreg", 4)),
  nm = var_names_miss
)

meth_miss

# Join with rest of variable names
meths <- setNames(rep("", ncol(dat_reg)), names(dat_reg))
meths[var_names_miss] <- meth_miss

meths


# Create matrix for MI
matpred <- matrix(1, ncol(dat_reg), ncol(dat_reg),
                  dimnames = list(names(dat_reg), names(dat_reg)))

# Don't impute a var using itself
diag(matpred) <- 0 

# Exclude outcomes
matpred[, c(
  "ci_s", 
  "ci_t", 
  "srv_t", 
  "srv_s", 
  "aa_auto" 
)] <- 0

matpred

# We don't impute any complete variable
matpred[!(rownames(matpred) %in% var_names_miss), ] <- 0

# Should read as: row gets imputed using cols with 1 as predictors
# If a row entirely 0: var is not going to be imputed
# If a col entirely 0: var is never used as predictor in imp of other vars
View(matpred)



# MICE --------------------------------------------------------------------


# Set for smcfcs too
m <- 10 # 10
iters <- 15 # 15

# Run imputations - NO RUN
imps_mice <- mice(
  data = dat_reg, 
  m = m, 
  maxit = iters,  
  method = meths, 
  predictorMatrix = matpred
)


plot(imps_mice)



# For relapse
mice_rel <- with(
  imps_mice,
  coxph(formula = Surv(ci_t, ci_s == 1) ~ classification_diag1 + age_allo1 + stage_cr_allo1 + 
          donor + source + cmvpat + ric + karnofsk_allo1 + cytog_recat + 
          tcelexvivo_allo1)
) %>% 
  pool() %>% 
  summary()


# For NRM
mice_nrm <- with(
  imps_mice,
  coxph(formula = Surv(ci_t, ci_s == 2) ~ classification_diag1 + age_allo1 + stage_cr_allo1 + 
          donor + source + cmvpat + karnofsk_allo1 + cytog_recat + 
          tcelexvivo_allo1)
) %>% 
  pool() %>% 
  summary()


# smcfcs ------------------------------------------------------------------



meths[which(meths == "polr")] <- "podds"


imps_smcfcs <- smcfcs::smcfcs(
  originaldata = dat_reg, 
  smtype = "compet", 
  smformula = c(
    "Surv(ci_t, ci_s == 1) ~ classification_diag1 + age_allo1 + stage_cr_allo1 + 
          donor + source + cmvpat + ric + karnofsk_allo1 + cytog_recat + 
          tcelexvivo_allo1",
    "Surv(ci_t, ci_s == 2) ~ classification_diag1 + age_allo1 + stage_cr_allo1 + 
          donor + source + cmvpat + karnofsk_allo1 + cytog_recat + 
          tcelexvivo_allo1"
  ),
  m = m, 
  numit = iters,
  method = meths
)


# Pool results
implist_smcfcs <- mitools::imputationList(imps_smcfcs$impDatasets)

# Relapse 
smcfcs_rel <- with(
  implist_smcfcs,
  coxph(formula = Surv(ci_t, ci_s == 1) ~ classification_diag1 + age_allo1 + stage_cr_allo1 + 
          donor + source + cmvpat + ric + karnofsk_allo1 + cytog_recat + 
          tcelexvivo_allo1)
) %>% 
  mitools::MIcombine() %>% 
  summary()


# NRM 
smcfcs_nrm <- with(
  implist_smcfcs,
  coxph(formula = Surv(ci_t, ci_s == 2) ~ classification_diag1 + age_allo1 + stage_cr_allo1 + 
          donor + source + cmvpat + karnofsk_allo1 + cytog_recat + 
          tcelexvivo_allo1)
) %>% 
  mitools::MIcombine() %>% 
  summary()



# Diagnostics of convergence
ests_df <- purrr::map_dfr(1:m, function(i) {
  as.data.frame(t(imps_smcfcs$smCoefIter[i, ,])) %>% 
    mutate(iter = 1:iters)
}, .id = "imp_dat") %>% 
  data.table()

# Set colnames for plotting 
cols_names <- c(
  "imp_dat", 
  paste0(rownames(smcfcs_rel), ".1"), 
  paste0(rownames(smcfcs_nrm), ".2"),
  "iter"
)

setnames(ests_df, new = cols_names)

melt.data.table(
  data = ests_df, 
  variable.name = "covar",
  value.name = "value", 
  id.vars = c("imp_dat","iter")
) %>% 
  ggplot(aes(iter, value, col = imp_dat)) +
  geom_line() +
  theme(legend.position = "none") +
  ggplot2::facet_wrap(~ covar) +
  theme_bw()




# Compare results ---------------------------------------------------------



# Relapse:
mod_rel
mice_rel
smcfcs_rel

names_preds <- names(mod_rel$coefficients)

res_rel <- cbind.data.frame(
  "CCA" = exp(round(mod_rel$coefficients, 2)),
  "mice" = exp(round(mice_rel$estimate, 2)),
  "smcfcs" = exp(round(smcfcs_rel$results, 2)),
  "CCA_se" = round(broom::tidy(mod_rel)$std.error, 2),
  "mice_se" = round(mice_rel$std.error, 2),
  "smcfcs_se" = round(smcfcs_rel$se, 2)
)

rownames(res_rel) <- names_preds
View(res_rel)


# NRM
mod_nrm
mice_nrm
smcfcs_nrm


res_nrm <- cbind.data.frame(
  "CCA" = exp(round(mod_nrm$coefficients, 2)),
  "mice" = exp(round(mice_nrm$estimate, 2)),
  "smcfcs" = exp(round(smcfcs_nrm$results, 2)),
  "CCA_se" = round(broom::tidy(mod_nrm)$std.error, 2),
  "mice_se" = round(mice_nrm$std.error, 2),
  "smcfcs_se" = round(smcfcs_nrm$se, 2)
)

rownames(res_nrm) <- names_preds
View(res_nrm)




# Function for plotting smcfcs convergence --------------------------------


forms <- c("Surv(ci_t, ci_s == 1) ~ classification_diag1 + age_allo1 + stage_cr_allo1 + 
          donor + source + cmvpat + ric + karnofsk_allo1 + cytog_recat + 
          tcelexvivo_allo1",
"Surv(ci_t, ci_s == 2) ~ classification_diag1 + age_allo1 + stage_cr_allo1 + 
          donor + source + cmvpat + karnofsk_allo1 + cytog_recat + 
          tcelexvivo_allo1")


# Will need to change if model specified for censoring
ggplot_smcfcs_converg <- function(imps_obj,
                                  smformula,
                                  dat) {
  
  # Extract meta data
  m <- dim(imps_obj$smCoefIter)[1]
  iters <- dim(imps_obj$smCoefIter)[3]
  
  # Number of comp risk models
  K <- length(smformula)
  
  # Get column names for smcoefiter
  coef_names_list <- lapply(X = 1:K, FUN = function(k) {
    rhs <- gsub(x = smformula[k], pattern = ".*~", replacement = "")
    
    model_mat <- model.matrix(
      object <- as.formula(paste0("~ 1 +", rhs)), 
      data <- dat
    )
    
    model_mat <- model_mat[, !(colnames(model_mat) %in% "(Intercept)")]
    
    coef_names_modk <- paste0(colnames(model_mat), ".", as.character(k))
  })
  
  # Unlist for names of smcoefiter
  coef_names <- unlist(coef_names_list)
  
  # Diagnostics of convergence\
  ests_list <- lapply(X = 1:m, function(i) {
    
    coef_dat <- as.data.frame(t(imps_obj$smCoefIter[i, ,]))
    coef_dat$iter <- 1:iters
    coef_dat$imp <- i 
    
    return(coef_dat)
  })
  
  ests_df <- data.table(do.call(rbind.data.frame, ests_list))
  
  # Set names
  setnames(ests_df, new = c(coef_names, "iter", "imp"))
  
  # Make plot
  p <- melt.data.table(
    data = ests_df, 
    variable.name = "covar",
    value.name = "value", 
    id.vars = c("imp","iter")
  ) %>% 
    ggplot(aes(iter, value, col = factor(imp))) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "none") +
    ggplot2::facet_wrap(~ covar) 
  
  return(p)
}


ggplot_smcfcs_converg(
  imps_obj = imps_smcfcs,
  smformula = forms,
  dat = dat_reg
)




