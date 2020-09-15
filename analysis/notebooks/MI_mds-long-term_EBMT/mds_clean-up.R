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
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # or here

knitr::opts_chunk$set(
  echo = FALSE, 
  out.width = "100%", 
  warning = FALSE, 
  message = FALSE
)

options(knitr.kable.NA = '', scipen=999)

devtools::load_all()

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
  gtsummary,
  survminer
)

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2)) 

# Reprod.
set.seed(1984)


# Reading in data ---------------------------------------------------------


# Read-in data with hctci
dat_hctci <- sjlabelled::read_spss(
  path = "analysis/data/raw_data/Final_MDS_Longterm_20200903_LK.sav" 
) %>% 
  
  # Clean variable names (i.e. bring lower case, with _ in between)
  janitor::clean_names() %>% 
  
  # Make into data.table
  setDT() %>% 
  
  # Subset aa_auto for matching and new vars
  .[, .(
    aa_auto,
    hctci,
    hctci_risk,
    agedonor_allo1_1
  )]


# Load-in cytogenetics data
load("analysis/data/raw_data/cytodata.RData")
data.table::setDT(cytodata) 
setnames(cytodata, new = c("aa_auto", "cyto_score", "vchromos"))

# Recat into three ordered levels
cytodata[, ':=' (
  vchromos = droplevels(vchromos),
  cytog_threecat = cut(
    cyto_score, 
    breaks = c(0, 3, 4, Inf),
    include.lowest = T, 
    ordered_result = T,
    labels = c("good(<=3)", "poor(4)", "very-poor(5)")
  )
)]

# Read-in original data from Schetelig and LdW (with published subset)
dat_orig <- sjlabelled::read_spss(
  path = "analysis/data/raw_data/MDS_longterm_R_180614.sav"
) %>% 
  
  # Clean variable names
  janitor::clean_names() %>%
  setDT() 


# Combine all three files via Reduce + merge 
dat_combined <- Reduce(
  f = function(x, y) merge.data.table(x, y, all.x = T, by = "aa_auto"),
  x = list(dat_orig, dat_hctci, cytodata)
)

# Clean env
rm(dat_hctci, cytodata)


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
  #"donsex_allo1_1", Only needed if match_allo1_1 does not work out
  "match_allo1_1",
  "age_allo1", 
  "mdsclass", 
  "donorrel", 
  "agedonor_allo1_1", 
  "karnofsk_allo1",
  "crnocr", 
  "cmv_combi_allo1_1",
  "cytog_threecat",
  "hctci_risk"
)


# Data / variable formatting ----------------------------------------------


# Factors that will use their SPSS value labels 
vars_val_labs <- c(
  "mdsclass",
  "crnocr",
  "patsex",
  "donorrel",
  "hctci_risk",
  "cmv_combi_allo1_1",
  #"donsex_allo1_1",
  "match_allo1_1"
)

# Make dataset ready for regression
dat_mds_reg <- data.table::copy(dat_combined) %>% 
  
  # Select variables, Remove one patient with wrong donor age (data error)
  .[aa_auto != 357439, ..vars_keep] %>% 
  
  # Set time in years, age allo and donor age in decades
  .[, ':=' (
    ci_allo1 = ci_allo1 / 12,
    srv_allo1 = srv_allo1 / 12,
    agedonor_allo1_decades = agedonor_allo1_1 / 10,
    age_allo1_decades = age_allo1 / 10
  )] %>% 
  
  # Admin censoring at 120 mo. (10 years) 
  .[, ':=' (
    ci_allo1 = ifelse(ci_allo1 >= 10, 10, ci_allo1),
    srv_allo1 = ifelse(srv_allo1 >= 10, 10, srv_allo1)
  )] %>% 
  
  # Adjust statuses based on admin cens
  .[, ':=' (
    ci_s_allo1 = ifelse(ci_allo1 == 10, 0, sjlabelled::as_numeric(ci_s_allo1)),
    srv_s_allo1 = ifelse(srv_allo1 == 10, 0, sjlabelled::as_numeric(srv_s_allo1))
  )] %>% 
  
  # For some vars, use the spss values labels as factor labels
  .[, (vars_val_labs) := lapply(.SD, sjlabelled::as_label), 
    .SDcols = vars_val_labs] %>% 
  
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
  
  # Reorder crnocr, relevel sex, order hctci
  .[, ':=' (
    crnocr = factor(
      crnocr, 
      levels = c("Untreated/not aimed at remission", "noCR", "CR")
    ),
    patsex = relevel(patsex, ref = "Female"),
    #donsex_allo1_1 = relevel(donsex_allo1_1, ref = "Female"),
    hctci_risk = as.ordered(hctci_risk)
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


# Table 1 -----------------------------------------------------------------

#' # Basic descriptives

dat_mds_reg %>% 
  .[, !c(
    "srv_s_allo1", 
    "srv_allo1",
    "aa_auto", 
    "ci_allo1"
  )] %>% 

  # Set label if desired
  tbl_summary(by = "mdsclass")# %>%
  #add_n(statistic = "{p_miss}%", col_label = "% missing")

# Visualise missingness
naniar::gg_miss_var(x = dat_mds_reg, facet = mdsclass, show_pct = T)
naniar::gg_miss_upset(dat_mds_reg, nsets = n_var_miss(dat_mds_reg))


# Univ cox models relapse -------------------------------------------------


#' # Univariable cox models relapse


#  Vector of predictors mod1 
predictors_mod1 <- c(
  "age_allo1_decades",
  "mdsclass", 
  "donorrel", 
  "agedonor_allo1_decades",
  "karnofsk_allo1", 
  "match_allo1_1",
  "crnocr", 
  "cmv_combi_allo1_1",
  "cytog_threecat",
  "hctci_risk"
)

names(predictors_mod1) <- predictors_mod1

# Univariate cause-specific models (relapse)
univ_mods_rel <- purrr::map(
  .x = predictors_mod1,
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
preds <- paste(predictors_mod1, collapse = " + ")
preds 

# For relapse
form_rel <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 1) ~ ", preds))
form_rel

mod_rel <- coxph(formula = form_rel, data = dat_mds_reg)


summary(mod_rel)
cox.zph(mod_rel)



# For nrm
form_nrm <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", preds))
form_nrm
mod_nrm <- coxph(formula = form_nrm, data = dat_mds_reg)

summary(mod_nrm)
cox.zph(mod_nrm)

#dev.new()
#ggforest(mod_rel, data = dat_mds_reg, main = "Rel")

#dev.new()
#ggforest(mod_nrm, data = dat_mds_reg, main = "NRM")


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


# Which vars actually have some data missing? (no binary vars)
var_names_miss <- naniar::miss_var_which(dat_mds_reg)

set_mi_methods <- function(dat,
                           var_names_miss,
                           imp_type = "mice",
                           cont_method = "norm") {
  
  # Set up methods vector
  n_vars_miss <- length(var_names_miss)
  meths_miss <- setNames(character(n_vars_miss), var_names_miss)
  
  # Get indicators:
  ordered_ind <- sapply(dat_mds_reg[, ..var_names_miss], is.ordered)
  contin_ind <- sapply(dat_mds_reg[, ..var_names_miss], is.numeric)
  
  unordered_ind <- sapply(
    dat_mds_reg[, ..var_names_miss], 
    function(col) length(levels(col)) > 2 & !is.ordered(col)
  )
  
  binary_ind <- sapply(
    dat_mds_reg[, ..var_names_miss], 
    function(col) length(levels(col)) == 2 & !is.ordered(col)
  )
  
  # Set up methods 
  meths_miss[ordered_ind] <- "polr"
  meths_miss[contin_ind] <- cont_method
  meths_miss[unordered_ind] <- "polyreg"
  meths_miss[binary_ind] <- "logreg"
  
  # Adjust if for smcfcs
  if (imp_type == "smcfcs") {
    meths_miss[which(meths_miss == "polyreg")] <- "mlogit"
    meths_miss[which(meths_miss == "polr")] <- "podds"
  }
  
  # Make global vec
  meths <- setNames(character(ncol(dat)), names(dat))
  meths[var_names_miss] <- meths_miss
  
  return(meths)
}


# Set methods accordingly
meths <- set_mi_methods(
  dat = dat_mds_reg,
  var_names_miss = var_names_miss, 
  imp_type = "mice",
  cont_method = "norm"
)


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
  "patsex"
)] <- 0


# We don't impute any complete variable
matpred[!(rownames(matpred) %in% var_names_miss), ] <- 0

# Should read as: row gets imputed using cols with 1 as predictors
# If a row entirely 0: var is not going to be imputed
# If a col entirely 0: var is never used as predictor in imp of other vars
#View(matpred)


# Run mice ----------------------------------------------------------------


# Set for smcfcs too
m <- 20
iters <- 15

# Run imputations 
imps_mice <- mice(
  data = dat_mds_reg,
  m = m,
  maxit = iters,
  method = meths,
  predictorMatrix = matpred
)

# Load object 
saveRDS(imps_mice, file = "imps_mice.rds")

# Diagnostic
plot(imps_mice)


# Combine results (will do this later in one go with mstate)

# For relapse
mice_rel <- with(
  imps_mice,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ age_allo1_decades + mdsclass + 
          donorrel + agedonor_allo1_decades + karnofsk_allo1 + match_allo1_1 + 
          crnocr + cmv_combi_allo1_1 + cytog_threecat + hctci_risk)
) %>% 
  pool() %>% 
  summary()

# For NRM
mice_nrm <- with(
  imps_mice,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ age_allo1_decades + mdsclass + 
          donorrel + agedonor_allo1_decades + karnofsk_allo1 + match_allo1_1 + 
          crnocr + cmv_combi_allo1_1 + cytog_threecat + hctci_risk)
) %>% 
  pool() %>% 
  summary()


# Run smcfcs --------------------------------------------------------------


# smcfcs will throw an error if there are variables in the data
# That will not be used at all (as per matpred), but still have missings

# Change coding of polr and polyreg, exclude the prev vars
meths <- set_mi_methods(
  dat = dat_mds_reg,
  var_names_miss = var_names_miss, 
  imp_type = "smcfcs",
  cont_method = "norm"
)

# No need for matpred smcfcs will use other var of the subs model 

# Make formula
smform_smcfcs <- c(
  Reduce(paste, deparse(form_rel)),
  Reduce(paste, deparse(form_nrm))
)

# You need the == ==
imps_smcfcs <- smcfcs::smcfcs(
  originaldata = dat_mds_reg,
  smtype = "compet",
  smformula = smform_smcfcs,
  m = m,
  numit = iters,
  method = meths
)

saveRDS(imps_smcfcs, file = "imps_smcfcs.rds")
#imps_smcfcs <- readRDS("imp_objects/imps_smcfcs.rds")

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
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ age_allo1_decades + mdsclass + 
          donorrel + agedonor_allo1_decades + karnofsk_allo1 + match_allo1_1 + 
          crnocr + cmv_combi_allo1_1 + cytog_threecat + hctci_risk)
) %>% 
  mitools::MIcombine() %>% 
  summary()


# NRM 
smcfcs_nrm <- with(
  implist_smcfcs,
  coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ age_allo1_decades + mdsclass + 
          donorrel + agedonor_allo1_decades + karnofsk_allo1 + match_allo1_1 + 
          crnocr + cmv_combi_allo1_1 + cytog_threecat + hctci_risk)
) %>% 
  mitools::MIcombine() %>% 
  summary()





# smcfcs with centre as aux -----------------------------------------------


# We want to use centre info as covar in imputation, but not in subst. model
# In smcfcs doing this is possible, but implies centre is indepdendent
# of outcome, conditional on covariates in the subst. model
# See: https://cran.r-project.org/web/packages/smcfcs/vignettes/smcfcs-vignette.html


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




# Forest attempts ---------------------------------------------------------


# Try relapse first
rel_smcfcs_tidy <- smcfcs_rel %>% 
  data.table(keep.rownames = "term") %>% 
  setnames(old = c("results", "se"), new = c("estimate", "std.error")) %>% 
  .[, .(term, estimate, std.error)]

dat_rel <- rbind(rel_smcfcs_tidy, mice_rel, broom::tidy(mod_rel),
                 fill = T, idcol = "method") %>% 
  .[, .(method, term, estimate, std.error)] %>% 
  .[, method := factor(
    method,
    levels = 1:3,
    labels = c("smcfcs", "mice", "CCA")
  )] %>% 
  
  # Make CI
  .[, ':=' (
    estimate = exp(estimate),
    CI_low = exp(estimate - pnorm(0.975) * std.error),
    CI_upp = exp(estimate + pnorm(0.975) * std.error)
  )] %>% 
  
  # colors
  .[, color_row := rep(c("white", "gray"), length.out = .N)]


p_rel <- dat_rel %>% 
  ggplot(aes(x = term, y = estimate, group = method, shape = method)) +
  geom_vline(aes(xintercept = term, col = color_row), size = 13) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_upp, xmin = term, xmax = term),
                 position = position_dodge(width = 0.75), size = 1) +
  geom_point(position = position_dodge2(width = 0.75), size = 2) +
  scale_y_continuous(trans = "log", breaks = c(0.5, 1, 2, 3, 5)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw(base_size = 14) + 
  coord_flip(ylim = c(0.5, 4)) +
  scale_color_manual(values = c("white", "lightgray"), guide = "none") + 
  ggtitle("Relapse") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.5)) + 
  #scale_color_discrete_qualitative(guide = "none", palette = "Dark 3")+
  xlab(NULL) +
  ylab("Hazard ratio (95% CI)")




# Try nrm 
nrm_smcfcs_tidy <- smcfcs_nrm %>% 
  data.table(keep.rownames = "term") %>% 
  setnames(old = c("results", "se"), new = c("estimate", "std.error")) %>% 
  .[, .(term, estimate, std.error)]

dat_nrm <- rbind(nrm_smcfcs_tidy, mice_nrm, broom::tidy(mod_nrm),
                 fill = T, idcol = "method") %>% 
  .[, .(method, term, estimate, std.error)] %>% 
  .[, method := factor(
    method,
    levels = 1:3,
    labels = c("smcfcs", "mice", "CCA")
  )] %>% 
  
  # Make CI
  .[, ':=' (
    estimate = exp(estimate),
    CI_low = exp(estimate - pnorm(0.975) * std.error),
    CI_upp = exp(estimate + pnorm(0.975) * std.error)
  )] %>% 
  
  # colors
  .[, color_row := rep(c("white", "gray"), length.out = .N)]


p_nrm <- dat_nrm %>% 
  ggplot(aes(x = term, y = estimate, group = method, shape = method)) +
  geom_vline(aes(xintercept = term, col = color_row), size = 13) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_upp, xmin = term, xmax = term),
                 position = position_dodge(width = 0.75), size = 1,
                 col = "black") +
  geom_point(position = position_dodge2(width = 0.75), size = 2) +
  scale_y_continuous(trans = "log", breaks = c(0.5, 1, 2, 3, 5)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw(base_size = 14) + 
  coord_flip(ylim = c(0.5, 4)) +
  scale_color_manual(values = c("white", "lightgray"), guide = "none") + 
  ggtitle("NRM") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(), 
        plot.title = element_text(hjust = 0.5)) + 
  xlab(NULL) +
  ylab("Hazard ratio (95% CI)") 

p_nrm


ggpubr::ggarrange(p_rel, p_nrm + theme(axis.text.y = element_blank()), 
                  common.legend = T, ncol = 2, widths = c(1.5, 1), 
                  legend = "bottom")
