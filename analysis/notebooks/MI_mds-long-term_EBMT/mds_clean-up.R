##***********************************##
## MDS clean-up (not incl. AFT mods) ##
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
  
  # aa_auto and covariates
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
  "tceldepl_allo1"
)

# Omitted possible covars:
# c("source_allo1", "tbi_allo1", "cmv_combi_allo1_1", "DONRL_allo1_1",
# "ric_allo1", "match_allo1_1", "centre_allo1)


# Data / variable formatting ----------------------------------------------


# Factors that will use their SPSS value labels 
vars_val_labs <- c(
  "mdsclass",
  "crnocr",
  "tceldepl_allo1",
  "patsex",
  "donorrel"
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
    agedonor_allo1_1 = agedonor_allo1_1 / 10,
    age_allo1 = age_allo1 / 10
  )] %>% 
  
  # For some vars, use the spss values labels as factor labels
  .[, (vars_val_labs) := lapply(.SD, sjlabelled::as_label), 
    .SDcols = vars_val_labs] %>% 
    
  # Group categories of t-cell depletion
  .[, tceldepl_allo1 := fcase(
    tceldepl_allo1 == "yes invivo, no exvivo", "invivo_only",
    tceldepl_allo1 == "yes invivo+exvivo", "exvivo",
    tceldepl_allo1 == "yes exvivo, no invivo", "exvivo",
    tceldepl_allo1 == "no", "no"
  )] %>% 
  
  # Make tcels into factor
  .[, tceldepl_allo1 := factor(
    tceldepl_allo1, 
    levels = c("no", "invivo_only", "exvivo")
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
  
  # Center year of transplant at zero, reorder crnocr
  .[, ':=' (
    crnocr = factor(
      crnocr, 
      levels = c("Untreated/not aimed at remission", "noCR", "CR")
    ),
    year_allo1 = year_allo1 - mean(year_allo1)
  )] %>% 
  
  # Dichtomise bm_allo1 at <= 1
  .[, dichot_bm_allo1 := cut(
    x = bm_allo1, 
    breaks = c(0, 1, 100),
    include.lowest = TRUE, 
    labels = c("<=1", ">1") 
  )]
  
# Drop unused factor levels 
which_factors <- names(dat_mds_reg)[sapply(dat_mds_reg, is.factor)]
dat_mds_reg[, (which_factors) := lapply(.SD, droplevels), .SDcols = which_factors]

# Conserve var labels from original data (for mutated vars)
dat_mds_reg <- data.table::copy(
  sjlabelled::copy_labels(df_new = dat_mds_reg, df_origin = dat_orig)
)

# View
dat_mds_reg


# Table 1 -----------------------------------------------------------------


dat_mds_reg %>% 
  .[, !c(
    "srv_s_allo1", 
    "srv_allo1",
    "aa_auto", 
    "pb_diag1", 
    "pb_allo1", 
    "bm_diag1",
    "ci_allo1"
  )] %>% 
  
  # Unorder karnofsky or it causes issues
  .[, karnofsk_allo1 := factor(karnofsk_allo1, ordered = F)] %>% 
  
  # Set label if desired
  tbl_summary(by = "mdsclass") %>%
  add_n(statistic = "{p_miss}%", col_label = "% missing")

# Visualise missingness
naniar::gg_miss_var(x = dat_mds_reg, facet = mdsclass, show_pct = T)
naniar::gg_miss_upset(dat_mds_reg, nsets = n_var_miss(dat_mds_reg))


# Density plots MDS / cross tables with CR --------------------------------


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


#  Vector of  predictors 
predictors <- c("patsex", "age_allo1",
                "mdsclass", "donorrel", 
                "dichot_bm_allo1", "agedonor_allo1_1",
                "karnofsk_allo1", "crnocr")

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


# RHS of formula
preds <- paste(predictors, collapse = " + ")
preds 

# For relapse
form_rel <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 1) ~ ", preds))
form_rel

mod_rel <- coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1 + mdsclass + 
                   donorrel + dichot_bm_allo1 + agedonor_allo1_1 + karnofsk_allo1 + 
                   crnocr, 
                 data = dat_mds_reg)


summary(mod_rel)
cox.zph(mod_rel)



# For nrm
form_nrm <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", preds))
form_nrm
mod_nrm <- coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1 + mdsclass + 
                   donorrel + dichot_bm_allo1 + agedonor_allo1_1 + karnofsk_allo1 + 
                   crnocr, 
                 data = dat_mds_reg)

summary(mod_nrm)
cox.zph(mod_nrm)


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
  object = c(rep("", 4), "norm", "polr", "polyreg", "", "logreg"),
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

# To exclude from all imputation models:
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
  "tceldepl_allo1"
)] <- 0

matpred

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
View(matpred)


# Run mice ----------------------------------------------------------------


# Set for smcfcs too
m <- 10
iters <- 15

# Run imputations
imps_mice <- mice(
  data = dat_mds_reg, 
  m = m, 
  maxit = iters,  
  method = meths, 
  predictorMatrix = matpred
)

# Diagnostic
plot(imps_mice)


# Combine results (will do this later in one go with mstate)

# For relapse
with(
  imps_mice,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1 + mdsclass + 
          donorrel + dichot_bm_allo1 + agedonor_allo1_1 + karnofsk_allo1 + 
          crnocr)
) %>% 
  pool() %>% 
  summary()

# Compare to CCA
mod_rel

# For NRM
with(
  imps_mice,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1 + mdsclass + 
          donorrel + dichot_bm_allo1 + agedonor_allo1_1 + karnofsk_allo1 + 
          crnocr)
) %>% 
  pool() %>% 
  summary()

# Compare to CCA
mod_nrm


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

# You need the == ==
imps_smcfcs <- smcfcs::smcfcs(
  originaldata = dat_mds_reg, 
  smtype = "compet", 
  smformula = c(
    "Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1 + mdsclass + 
    donorrel + dichot_bm_allo1 + agedonor_allo1_1 + karnofsk_allo1 + 
    crnocr",
    "Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1 + mdsclass + 
    donorrel + dichot_bm_allo1 + agedonor_allo1_1 + karnofsk_allo1 + 
    crnocr"
  ),
  m = 10, 
  numit = 15,
  method = meths
)


# 
# Pool results
implist_smcfcs <- mitools::imputationList(imps_smcfcs$impDatasets)

# Relapse 
rel_smcfcs <- with(
  implist_smcfcs,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ patsex + age_allo1 + mdsclass + 
          donorrel + dichot_bm_allo1 + agedonor_allo1_1 + karnofsk_allo1 + 
          crnocr)
) %>% 
  mitools::MIcombine() %>% 
  summary()

mod_rel

# NRM 
with(
  implist_smcfcs,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ patsex + age_allo1 + mdsclass + 
          donorrel + dichot_bm_allo1 + agedonor_allo1_1 + karnofsk_allo1 + 
          crnocr)
) %>% 
  mitools::MIcombine() %>% 
  summary()

mod_nrm

# Diagnostics of convergence
ests_df <- purrr::map_dfr(1:m, function(i) {
  as.data.frame(t(imps_smcfcs$smCoefIter[i, ,])) %>% 
    mutate(iter = 1:iters)
}, .id = "imp_dat") %>% 
  data.table()

# Set colnames for plotting 
cols_names <- c(
  "imp_dat", 
  paste0(rownames(rel_smcfcs), ".1"), 
  paste0(rownames(rel_smcfcs), ".2"),
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
  ggplot2::facet_wrap(~ covar)


