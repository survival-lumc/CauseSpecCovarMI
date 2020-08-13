#' ---
#' title: "MDS analysis"
#' author: "Ed Bonneville"
#' date: "`r Sys.setenv(LANG = 'en_US.UTF-8'); format(Sys.Date(), '%d %B %Y')`"
#' output:
#'   html_document:
#'     df_print: kable
#'     toc: yes
#'     toc_float:
#'       collapsed: no
#'       smooth_scroll: no
#' always_allow_html: yes
#' ---


#+ echo = FALSE, message = FALSE, warning = FALSE
# SETUP -----------------------------------------------------------------------
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(
  echo = FALSE, 
  out.width = "100%", 
  warning = FALSE, 
  message = FALSE
)
options(knitr.kable.NA = '', scipen=999)



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
  broom,
  naniar,
  ggpubr, # will load ggplot with it
  gtsummary
)

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2)) 

# Reprod.
set.seed(1984)


# Reading in data ---------------------------------------------------------


# Read-in data with extra vars
dat_blasts <- sjlabelled::read_spss(
  path = "analysis/data/raw_data/Final_MDS_Longterm_20200106_LK.sav" 
) %>% 
  
  # Clean variable names (i.e. bring lower case, with _ in between)
  janitor::clean_names() %>% 
  
  # Make into data.table
  data.table() %>% 
  
  # Subset aa_auto for matching and new vars
  .[, .(
    aa_auto,
    bm_diag1,
    bm_allo1,
    pb_diag1,
    pb_allo1,
    agedonor_allo1_1
  )]


# Read-in original data for Schetelig and LdW (with published subset)
dat_orig <- sjlabelled::read_spss(
  path = "analysis/data/raw_data/MDS_longterm_R_180614.sav"
) %>% 
  
  # Clean variable names
  janitor::clean_names() %>%
  data.table() %>% 
  
  # Merge data - warning is just because of 'display_width' attr - can ignore
  data.table::merge.data.table(
    x = ., 
    y = dat_blasts, 
    all.x = TRUE, 
    by = "aa_auto"
  )


# Clean env
rm(dat_blasts)


# Define variables to keep - see table for motivation
vars_keep <- c(
  
  # Outcome vars
  "ci_s_allo1", 
  "ci_allo1", 
  "srv_s_allo1", 
  "srv_allo1",
  
  # aa_auto and (possible) covariates
  "aa_auto", 
  "patsex", 
  "age_allo1", 
  "year_allo1",
  "mdsclass", 
  "donorrel", 
  "bm_allo1", 
  "pb_allo1", # take it ref. DJE paper
  "pb_diag1", 
  "bm_diag1", 
  "agedonor_allo1_1", 
  "karnofsk_allo1",
  "crnocr", 
  "tceldepl_allo1",
  "centre_allo1",
  "vcmvpat_allo1"
)

# Omitted possible covars:
# c("source_allo1", "tbi_allo1", "DONRL_allo1_1",
# "ric_allo1", "match_allo1_1")


# Data / variable formatting ----------------------------------------------


# Factors that will use their SPSS value labels 
vars_val_labs <- c(
  "mdsclass",
  "crnocr",
  "tceldepl_allo1",
  "patsex",
  "donorrel",
  "centre_allo1",
  "vcmvpat_allo1"
)

# Make dataset ready for regression
dat_mds_reg <- data.table::copy(dat_orig) %>% 
  
  # Select variables
  .[, ..vars_keep] %>% 
  
  # Set time in years
  .[, ':=' (
    ci_allo1 = ci_allo1 / 12,
    srv_allo1 = srv_allo1 / 12
  )] %>% 
  
  # Admin censoring at 120 mo., age allo and donor age in decades
  .[, ':=' (
    ci_allo1 = ifelse(ci_allo1 >= 10, 10, ci_allo1),
    ci_s_allo1 = ifelse(ci_allo1 == 10, 0, sjlabelled::as_numeric(ci_s_allo1)),
    agedonor_allo1_decades = agedonor_allo1_1 / 10,
    age_allo1_decades = age_allo1 / 10
  )] %>% 
  
  # For some vars, use the spss values labels as factor labels
  .[, (vars_val_labs) := lapply(.SD, sjlabelled::as_label), 
    .SDcols = vars_val_labs] %>% 
  
  # cmv patients missings code
  .[vcmvpat_allo1 %in% c("unknown", "Not evaluated"), 
    vcmvpat_allo1 := NA_character_] %>% 
    
  # Group categories of t-cell depletion
  .[, tceldepl_allo1 := fcase(
    tceldepl_allo1 == "yes invivo, no exvivo", "invivo_only",
    tceldepl_allo1 == "yes invivo+exvivo", "exvivo",
    tceldepl_allo1 == "yes exvivo, no invivo", "exvivo",
    tceldepl_allo1 == "no", "no"
  )] %>% 
  
  # Make tcels into ordered factor
  .[, tceldepl_allo1 := factor(
    tceldepl_allo1, 
    levels = c("no", "invivo_only", "exvivo"), 
    ordered = T
  )] %>% 
  
  # Make Karnofsky into ordered three cats
  .[, karnofsk_allo1 := cut(
    x = karnofsk_allo1, 
    breaks = c(-Inf, 70, 80, 100), 
    include.lowest = T, 
    ordered_result = T, 
    labels = c("<=70", "80", ">=90")
  )] %>% 
  
  # Reverse ordering to have >= 90 as reference
  .[, karnofsk_allo1 := factor(
    karnofsk_allo1, levels = rev(levels(karnofsk_allo1))
  )] %>% 
  
  # Reorder crnocr, relevel sex
  .[, ':=' (
    crnocr = factor(
      crnocr, 
      levels = c("Untreated/not aimed at remission", "noCR", "CR")
    ),
    patsex = relevel(patsex, ref = "Female")
  )] %>% 
  
  # Dichtomise bm_allo1 at <= 1
  .[, dichot_pb_allo1 := cut(
    x = pb_allo1, 
    breaks = c(0, 1, 100),
    include.lowest = TRUE, 
    labels = c("<=1", ">1") 
  )] %>% 
  
  # Delete redundant age and donor age not in decades (they are in dat_orig)
  .[, ':=' (
    agedonor_allo1_1 = NULL,
    age_allo1 = NULL
  )]
  
# Drop unused factor levels 
which_factors <- names(dat_mds_reg)[sapply(dat_mds_reg, is.factor)]
dat_mds_reg[, (which_factors) := lapply(.SD, droplevels), .SDcols = which_factors]

# Conserve var labels from original data (for mutated vars)
dat_mds_reg <- data.table::copy(
  sjlabelled::copy_labels(df_new = dat_mds_reg, df_origin = dat_orig)
)

# View
#dat_mds_reg


# Table 1 -----------------------------------------------------------------

#' # Basic descriptives

dat_mds_reg %>% 
  .[, !c(
    "srv_s_allo1", 
    "srv_allo1",
    "aa_auto", 
    "pb_diag1", 
    #"pb_allo1", 
    "bm_diag1",
    "ci_allo1",
    "centre_allo1"
  )] %>% 
  
  # Unorder karnofsky and tcels or it causes issues
  .[, ':=' (
    karnofsk_allo1 = factor(karnofsk_allo1, ordered = F),
    tceldepl_allo1 = factor(tceldepl_allo1, ordered = F)
  )] %>% 

  # Set label if desired
  tbl_summary(by = "mdsclass") %>%
  add_n(statistic = "{p_miss}%", col_label = "% missing")

# Visualise missingness
naniar::gg_miss_var(x = dat_mds_reg, facet = mdsclass, show_pct = T)
naniar::gg_miss_upset(dat_mds_reg, nsets = n_var_miss(dat_mds_reg))


# Density plots MDS / cross tables with CR --------------------------------

#' # Density plots MDS / cross tables with CR


# Select blast variables
blast_vars <- stringr::str_subset(names(dat_mds_reg), "^bm_|^pb_")

list_densplots <- purrr::map(
  .x = blast_vars,
  .f = ~ {
    dat_mds_reg %>% 
      ggplot(aes(x = .[[.x]], fill = mdsclass)) +
      geom_density(alpha= .5, adjust = 1, na.rm = T) +
      xlim(c(0, 50)) +
      ggtitle(.x) +
      theme(axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) 
  }
)

p <- ggarrange(plotlist = list_densplots, nrow = 2, ncol = 2, 
               common.legend = T, legend = "top")

annotate_figure(p, left = "Density", bottom = "Blast percentage")


# For CR groups
blast_vars_allo <- stringr::str_subset(names(dat_mds_reg), "^bm_allo1|^pb_allo1")

list_densplots_allo <- purrr::map(
  .x = blast_vars_allo,
  .f = ~ {
    dat_mds_reg %>% 
      .[!is.na(crnocr)] %>% 
      ggplot(aes(x = .[[.x]], fill = crnocr)) +
      geom_density(alpha= .5, adjust = 1, na.rm = T) +
      xlim(c(0, 50)) +
      ggtitle(.x) +
      theme(axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      scale_colour_discrete(na.translate = F)
  }
)

p_allo <- ggarrange(
  plotlist = list_densplots_allo,
  ncol = 2, 
  common.legend = T, 
  legend = "top"
)

annotate_figure(p_allo, left = "Density", bottom = "Blast percentage (at allo)")


# Cross tables crnocr and mds class
dat_mds_reg %>% 
  .[, c("crnocr", "mdsclass")] %>% 
  tbl_summary(by = mdsclass)

dat_mds_reg %>% 
  .[, c("crnocr", "bm_allo1")] %>% 
  tbl_summary(by = "crnocr")


# Univ cox models relapse -------------------------------------------------


#' # Univariable cox models relapse



#  Vector of  predictors 
predictors <- c("patsex", "age_allo1_decades",
                "mdsclass", "donorrel", 
                "dichot_pb_allo1", "agedonor_allo1_decades",
                "karnofsk_allo1", "crnocr", "vcmvpat_allo1")

names(predictors) <- predictors

# Univariate cause-specific models (relapse)
univ_mods_rel <- purrr::map(
  .x = predictors,
  .f = ~ {
    form <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 1) ~ ", .x))
    mod <- coxph(form, data = dat_mds_reg)
  }
)

# Print any summary
summary(univ_mods_rel$karnofsk_allo1)
plot(cox.zph(univ_mods_rel$karnofsk_allo1))

# Or all
purrr::map(univ_mods_rel, summary)

# Plot all coxzphs
purrr::map(univ_mods_rel, ~ plot(cox.zph(.x)))



# Univ cox models nrm -----------------------------------------------------

#' # Univariable cox models nrm


# Univariate cause-specific models (nrm)
univ_mods_nrm <- purrr::map(
  .x = predictors,
  .f = ~ {
    form <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", .x))
    mod <- coxph(form, data = dat_mds_reg)
  }
)

# Print any summary
summary(univ_mods_nrm$karnofsk_allo1)
plot(cox.zph(univ_mods_nrm$karnofsk_allo1))

# Or all
purrr::map(univ_mods_nrm, summary)

# Plot all coxzphs
purrr::map(univ_mods_nrm, ~ plot(cox.zph(.x)))




# Multivariable cox models ------------------------------------------------


#' # Multivariable cox models (CCA)


# RHS of formula
preds <- paste(predictors, collapse = " + ")
preds 

# For relapse
form_rel <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 1) ~ ", preds))
form_rel

mod_rel <- coxph(formula = form_rel, 
                 data = dat_mds_reg)


summary(mod_rel)
cox.zph(mod_rel)



# For nrm
form_nrm <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", preds))
form_nrm
mod_nrm <- coxph(formula = form_nrm, 
                 data = dat_mds_reg)

summary(mod_nrm)
cox.zph(mod_nrm)



# Multivariable models without pb_allo1 dichot ----------------------------



mod_rel_nopb <- update(mod_rel, . ~ . - dichot_pb_allo1)
mod_nrm_nopb <- update(mod_nrm, . ~ . - dichot_pb_allo1)


# Prep mice matrix --------------------------------------------------------


# First, for competing risks, we make indicators for comp events and
# compute cumulative hazards

# Get event indictors
dat_mds_reg[, ':=' (
  ev1 = ifelse(ci_s_allo1 == 1, 1, 0),
  ev2 = ifelse(ci_s_allo1 == 2, 1, 0)
)]

# Compute cumulative hazards
dat_mds_reg[, ':=' (
  H1 = mice::nelsonaalen(data = dat_mds_reg, timevar = ci_allo1, statusvar = ev1),
  H2 = mice::nelsonaalen(data = dat_mds_reg, timevar = ci_allo1, statusvar = ev2)
)]


# Which vars actually have some data missing?
var_names_miss <- naniar::miss_var_which(dat_mds_reg)

# Set methods accordingly
meth_miss <- setNames(
  object = c(rep("", 4),  "polr", "polyreg", "", "logreg", "norm", "logreg"),
  nm = var_names_miss
)

meth_miss

# Join with rest of variable names
meths <- setNames(rep("", ncol(dat_mds_reg)), names(dat_mds_reg))
meths[var_names_miss] <- meth_miss

meths


# Create matrix for MI
matpred <- matrix(1, ncol(dat_mds_reg), ncol(dat_mds_reg),
                  dimnames = list(names(dat_mds_reg), names(dat_mds_reg)))

# Don't impute a var using itself
diag(matpred) <- 0 

# To exclude from all imputation models as predictors:
matpred[, c(
  "ci_s_allo1", 
  "ci_allo1", 
  "srv_allo1", 
  "srv_s_allo1", 
  "aa_auto", 
  "year_allo1",
  "bm_allo1",
  "bm_diag1",
  "pb_allo1",
  "pb_diag1",
  "tceldepl_allo1",
  "centre_allo1" # exclude centres here
)] <- 0


# We don't impute any complete variable
matpred[!(rownames(matpred) %in% var_names_miss), ] <- 0

# We also don't impute the bm and pb (non dichot vars, along with tcel)
matpred[rownames(matpred) %in% c(
  "bm_allo1",
  "bm_diag1",
  "pb_allo1",
  "pb_diag1",
  "tceldepl_allo1"
), ] <- 0

# Should read as: row gets imputed using cols with 1 as predictors
# If a row entirely 0: var is not going to be imputed
# If a col entirely 0: var is never used as predictor in imp of other vars
#View(matpred)


# Run mice ----------------------------------------------------------------


# Set for smcfcs too
m <- 10 
iters <- 15 

# Run imputations 
# imps_mice <- mice(
#   data = dat_mds_reg, 
#   m = m, 
#   maxit = iters,  
#   method = meths, 
#   predictorMatrix = matpred
# )

# Load object 
imps_mice <- readRDS("imps_mice.rds")

# Diagnostic
plot(imps_mice)


# Combine results (will do this later in one go with mstate)

# For relapse
mice_rel <- with(
  imps_mice,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1_decades + mdsclass + 
          donorrel + dichot_pb_allo1 + agedonor_allo1_decades + karnofsk_allo1 + 
          crnocr + vcmvpat_allo1)
) %>% 
  pool() %>% 
  summary()

# For NRM
mice_nrm <- with(
  imps_mice,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1_decades + mdsclass + 
          donorrel + dichot_pb_allo1 + agedonor_allo1_decades + karnofsk_allo1 + 
          crnocr + vcmvpat_allo1)
) %>% 
  pool() %>% 
  summary()



# Center as aux var (with dichot_pb_allo1 in model) -----------------------


matpred[rowSums(matpred) > 0, "centre_allo1"] <- 1

# This will be much, much slower 
# imps_mice_centre <- mice(
#   data = dat_mds_reg, 
#   m = m, 
#   maxit = iters,  
#   method = meths, 
#   predictorMatrix = matpred
# )

# For relapse
# mice_rel_centre <- with(
#   imps_mice_centre,
#   coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1_decades + mdsclass + 
#           donorrel + dichot_pb_allo1 + agedonor_allo1_decades + karnofsk_allo1 + 
#           crnocr + vcmvpat_allo1)
# ) %>% 
#   pool() %>% 
#   summary()
# 
# # For NRM
# mice_nrm_centre <- with(
#   imps_mice_centre,
#   coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1_decades + mdsclass + 
#           donorrel + dichot_pb_allo1 + agedonor_allo1_decades + karnofsk_allo1 + 
#           crnocr + vcmvpat_allo1)
# ) %>% 
#   pool() %>% 
#   summary()


# MICE (no dichot_pb_allo1) -----------------------------------------------


# Exclude it from imps, exclude center again
matpred["dichot_pb_allo1", ] <- 0
matpred[, "dichot_pb_allo1"] <- 0
matpred[, "centre_allo1"] <- 0


# imps_mice_nopb <- mice(
#   data = dat_mds_reg, 
#   m = m, 
#   maxit = iters,  
#   method = meths, 
#   predictorMatrix = matpred
# )

# Load object 
imps_micenopb <- readRDS("imps_mice.rds")


# For relapse
mice_rel_nopb <- with(
  imps_mice,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1_decades + mdsclass + 
          donorrel  + agedonor_allo1_decades + karnofsk_allo1 + 
          crnocr + vcmvpat_allo1)
) %>% 
  pool() %>% 
  summary()

# For NRM
mice_nrm_nopb <- with(
  imps_mice,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1_decades + mdsclass + 
          donorrel + agedonor_allo1_decades + karnofsk_allo1 + 
          crnocr + vcmvpat_allo1)
) %>% 
  pool() %>% 
  summary()


# smcfcs: function to plot convergence ------------------------------------

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


# Run smcfcs --------------------------------------------------------------


# smcfcs will throw an error if there are variables in the data
# That will not be used at all (as per matpred), but still have missings
# (Like with bm_allo1, pb_allo1 etc.), let's exclude them
vars_excl <- c("bm_allo1", "pb_allo1", "bm_diag1", "pb_diag1", "tceldepl_allo1")
dat_mds_reg[, (vars_excl) := NULL]

# Change coding of polr and polyreg, exclude the prev vars
meths[which(meths == "polyreg")] <- "mlogit"
meths[which(meths == "polr")] <- "podds"
meths <- meths[!(names(meths) %in% vars_excl)]

# No need for matpred smcfcs will use other var of the subs model 

# Make formula
smform_smcfcs <- c(
  "Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1_decades + mdsclass +
    donorrel + dichot_pb_allo1 + agedonor_allo1_decades + karnofsk_allo1 +
    crnocr + vcmvpat_allo1",
  "Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1_decades + mdsclass +
    donorrel + dichot_pb_allo1 + agedonor_allo1_decades + karnofsk_allo1 +
    crnocr + vcmvpat_allo1"
)

# You need the == ==
# imps_smcfcs <- smcfcs::smcfcs(
#   originaldata = dat_mds_reg, 
#   smtype = "compet", 
#   smformula = smform_smcfcs,
#   m = m, 
#   numit = iters,
#   method = meths
# )

#saveRDS(imps_smcfcs, file = "imps_smcfcs.rds")
imps_smcfcs <- readRDS("imps_smcfcs.rds")

# Check convergence
ggplot_smcfcs_converg(
  imps_obj = imps_smcfcs, 
  smformula = smform_smcfcs, 
  dat = dat_mds_reg
)


# Pool results
implist_smcfcs <- mitools::imputationList(imps_smcfcs$impDatasets)

# Relapse 
smcfcs_rel <- with(
  implist_smcfcs,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1_decades + mdsclass + 
          donorrel + dichot_pb_allo1 + agedonor_allo1_decades + karnofsk_allo1 + 
          crnocr + vcmvpat_allo1)
) %>% 
  mitools::MIcombine() %>% 
  summary()


# NRM 
smcfcs_nrm <- with(
  implist_smcfcs,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1_decades + mdsclass + 
          donorrel + dichot_pb_allo1 + agedonor_allo1_decades + karnofsk_allo1 + 
          crnocr + vcmvpat_allo1)
) %>% 
  mitools::MIcombine() %>% 
  summary()





# smcfcs with centre as aux -----------------------------------------------


# We want to use centre info as covar in imputation, but not in subst. model
# In smcfcs doing this is possible, but implies centre is indepdendent
# of outcome, conditional on covariates in the subst. model
# See: https://cran.r-project.org/web/packages/smcfcs/vignettes/smcfcs-vignette.html


# # Exclude from df
# matpred_smcfcs <- matrix(0, ncol(dat_mds_reg), ncol(dat_mds_reg),
#                   dimnames = list(names(dat_mds_reg), names(dat_mds_reg)))
# 
# # Use covariates from subst as predictors
# matpred_smcfcs[, colnames(matpred_smcfcs) %in% predictors] <- 1
# diag(matpred_smcfcs) <- 0
# 
# # Add centre as aux
# matpred_smcfcs[, "centre_allo1"] <- 1
# 
# # We don't impute any complete variable
# matpred_smcfcs[!(rownames(matpred_smcfcs) %in% var_names_miss), ] <- 0
# 
# # Run
# # imps_smcfcs_centre <- smcfcs::smcfcs(
# #   originaldata = dat_mds_reg, 
# #   smtype = "compet", 
# #   smformula = smform_smcfcs,
# #   m = m, 
# #   numit = iters,
# #   method = meths,
# #   predictorMatrix = matpred_smcfcs
# # )
# 
# 
# 
# 
# # Pool results
# implist_smcfcs_centre <- mitools::imputationList(imps_smcfcs_nopb$impDatasets)
# 
# # Relapse 
# smcfcs_rel_centre <- with(
#   implist_smcfcs_centre,
#   coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1_decades + mdsclass + 
#           donorrel + agedonor_allo1_decades + karnofsk_allo1 + 
#           crnocr + vcmvpat_allo1)
# ) %>% 
#   mitools::MIcombine() %>% 
#   summary()
# 
# 
# # NRM 
# smcfcs_nrm_centre <- with(
#   implist_smcfcs_centre,
#   coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1_decades + mdsclass + 
#           donorrel + agedonor_allo1_decades + karnofsk_allo1 + 
#           crnocr + vcmvpat_allo1)
# ) %>% 
#   mitools::MIcombine() %>% 
#   summary()



# smcfcs no pb ------------------------------------------------------------

# Exclude from df
dat_mds_reg[, dichot_pb_allo1 := NULL]
meths <- meths[!(names(meths) %in% "dichot_pb_allo1")]


# Make formula without pb allo
smform_smcfcs_nopb <- c(
  "Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1_decades + mdsclass +
    donorrel + agedonor_allo1_decades + karnofsk_allo1 +
    crnocr + vcmvpat_allo1",
  "Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1_decades + mdsclass +
    donorrel + agedonor_allo1_decades + karnofsk_allo1 +
    crnocr + vcmvpat_allo1"
)

# Omit pb_allo from subs model
# imps_smcfcs_nopb <- smcfcs::smcfcs(
#   originaldata = dat_mds_reg, 
#   smtype = "compet", 
#   smformula = smform_smcfcs_nopb,
#   m = m, 
#   numit = iters,
#   method = meths
# )

imps_smcfcs_nopb <- readRDS("imps_smcfcs_nopb.rds")

# Check convergence
ggplot_smcfcs_converg(
  imps_obj = imps_smcfcs_nopb, 
  smformula = smform_smcfcs_nopb, 
  dat = dat_mds_reg
)

# Pool results
implist_smcfcs_nopb <- mitools::imputationList(imps_smcfcs_nopb$impDatasets)

# Relapse 
smcfcs_rel_nopb <- with(
  implist_smcfcs_nopb,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1_decades + mdsclass + 
          donorrel + agedonor_allo1_decades + karnofsk_allo1 + 
          crnocr + vcmvpat_allo1)
) %>% 
  mitools::MIcombine() %>% 
  summary()


# NRM 
smcfcs_nrm_nopb <- with(
  implist_smcfcs_nopb,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1_decades + mdsclass + 
          donorrel + agedonor_allo1_decades + karnofsk_allo1 + 
          crnocr + vcmvpat_allo1)
) %>% 
  mitools::MIcombine() %>% 
  summary()





# Comparison models with dichot_pb_allo1 ----------------------------------


#' # Comparison models (with dichot pb allo1)


# Relapse:
#mod_rel
#mice_rel
#smcfcs_rel

names_preds <- names(mod_rel$coefficients)

res_rel <- cbind.data.frame(
  "CCA" = mod_rel$coefficients,
  "mice" = mice_rel$estimate,
  "smcfcs" = smcfcs_rel$results,
  "CCA_se" = broom::tidy(mod_rel)$std.error,
  "mice_se" = mice_rel$std.error,
  "smcfcs_se" = smcfcs_rel$se
)

rownames(res_rel) <- names_preds
round(res_rel, 2) %>% 
  as.data.frame() %>% 
  knitr::kable(caption = "REL results")


# NRM
# mod_nrm
# mice_nrm
# smcfcs_nrm


res_nrm <- cbind.data.frame(
  "CCA" = mod_nrm$coefficients,
  "mice" = mice_nrm$estimate,
  "smcfcs" = smcfcs_nrm$results,
  "CCA_se" = broom::tidy(mod_nrm)$std.error,
  "mice_se" = mice_nrm$std.error,
  "smcfcs_se" = smcfcs_nrm$se
)

rownames(res_nrm) <- names_preds
round(res_nrm, 2) %>% 
  as.data.frame() %>% 
  knitr::kable(caption = "NRM results")




# Comparison models without dichot_pb_allo1 -------------------------------

#' # Comparison models (no dichot pb allo1)


names_preds <- names(mod_rel_nopb$coefficients)

res_rel_nopb <- cbind.data.frame(
  "CCA" = mod_rel_nopb$coefficients,
  "mice" = mice_rel_nopb$estimate,
  "smcfcs" = smcfcs_rel_nopb$results,
  "CCA_se" = broom::tidy(mod_rel_nopb)$std.error,
  "mice_se" = mice_rel_nopb$std.error,
  "smcfcs_se" = smcfcs_rel_nopb$se
)

rownames(res_rel_nopb) <- names_preds
round(res_rel_nopb, 2) %>% 
  as.data.frame() %>% 
  knitr::kable(caption = "REL results")


# NRM
res_nrm_pb <- cbind.data.frame(
  "CCA" = mod_nrm_nopb$coefficients,
  "mice" = mice_nrm_nopb$estimate,
  "smcfcs" = smcfcs_nrm_nopb$results,
  "CCA_se" = broom::tidy(mod_nrm_nopb)$std.error,
  "mice_se" = mice_nrm_nopb$std.error,
  "smcfcs_se" = smcfcs_nrm_nopb$se
)

rownames(res_nrm_pb) <- names_preds
round(res_nrm_pb, 2) %>% 
  as.data.frame() %>% 
  knitr::kable(caption = "NRM results")


