##***********************************##
## MDS clean-up (not incl. AFT mods) ##
##***********************************##

# Load libraries
devtools::load_all()

pacman::p_load(
  foreign,
  tidyverse,
  survival, 
  mstate,
  mice, 
  smcfcs,
  mitools,
  naniar,
  ggpubr,
  gtsummary
)

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2)) 

# Reprod.
set.seed(1984)


# Reading in data ---------------------------------------------------------


# Original dataset
dat_orig <- read.spss(
  file = "analysis/data/raw_data/Final_MDS_Longterm_20200106_LK.sav",
  to.data.frame = T,
  use.value.labels = T,
  use.missings = T, 
  max.value.labels = 8
)

# To run liesbeth script
mdsdata <-spss.get("analysis/data/raw_data/MDS_longterm_R_180614.sav", 
                   datevars=c('datallo1', 'datrel_1_allo1', 
                              'datallo2', 'datlast'), lowernames=TRUE,
                   max.value.labels = 25)

# Subset of pats in paper
LdW_schetelig_subset <- read.spss(
  file = "analysis/data/raw_data/mds_longterm_aaauto.sav",
  to.data.frame = T
)

# Centre data
centre_dat <- read.spss(
  file = "analysis/data/raw_data/MDS_longterm_R_180614.sav",
  to.data.frame = T,
  use.value.labels = T,
  use.missings = T
) %>% 
  select_keeplabs("AA_AUTO", "CENTRE_allo1")

# Variables to keep
vars_keep <- c("AA_AUTO", 
               "PATSEX", 
               "match_allo1_1", # instead of pat sex
               "age_allo1", "YEAR_allo1",
               "mdsclass", "donorrel", 
               #"ric_allo1", # start without
               "bm_allo1", # try one of the blast vars
               "pb_allo1", # take it ref. DJE paper
               "pb_diag1", "bm_diag1", 
               "AGEDONOR_allo1_1", 
               #"rel_s_allo1", "rel_allo1", 
               "ci_s_allo1", "ci_allo1", 
               "srv_s_allo1", "srv_allo1", #"DONRL_allo1_1", 
               "KARNOFSK_allo1",
               #"cmv_combi_allo1_1", 
               #tbi_allo1
               "crnocr", "tceldepl_allo1")
               #"source_allo1") # more optional

# Check code for missing relapse of 80 patients or so (check liesbeth code)

# Others: 
# crnocr
# tceldepl_allo1 

dat

# Keep variables from study 
dat_mds <- dat_orig[dat_orig$AA_AUTO %in% LdW_schetelig_subset$AA_AUTO,]

dat_mds <- dat_orig %>% 
  select_keeplabs(vars_keep) %>% 
  filter(AA_AUTO %in% LdW_schetelig_subset$AA_AUTO) %>% 
  left_join(centre_dat)

mdsdata$AA_AUTO <- mdsdata$aa.auto 

library(data.table)
# Check missings
dat_mds %>% 
  #setDT() %>% 
  #.[is.na(ci_allo1), ]
  filter(is.na(ci_allo1)) %>% 
  select(AA_AUTO, ci_allo1, ci_s_allo1, srv_allo1) %>% 
  left_join(mdsdata, by = "AA_AUTO") %>% 
  select(AA_AUTO, ci.allo1, srv.allo1, ci.s.allo1, srv_allo1, ci_allo1) %>% 
  arrange(AA_AUTO)


dat_mds %>% 
  group_by(source_allo1) %>% 
  summarize(age = median(age_allo1))

dat_mds %>% 
  filter(AGEDONOR_allo1_1 < 18) %>% 
  View()

dat_mds$source_allo1
# Prepare rest of data ----------------------------------------------------


naniar::gg_miss_var(
  x = dat_mds,
  show_pct = T
)

dat_prepped <- dat_mds %>% 
  
  # In years or not? (currently in months)
  mutate(time = ci_allo1) %>% 
  
  # Exclude 86 patients with missing relapse info
  filter(!is.na(time)) %>% 
  
  # Create indicators, admin cens at 10 years
  mutate(
    delta = ifelse(time == 10 * 12, 0, as.numeric(ci_s_allo1) - 1), # should be factor in imp
    ev1 = ifelse(delta == 1, 1, 0),
    ev2 = ifelse(delta == 2, 1, 0)
  ) %>% 
  
  # Prep cum hazards
  mutate(
    H1 = nelsaalen_timefixed(., timevar = time, statusvar = ev1),
    H2 = nelsaalen_timefixed(., timevar = time, statusvar = ev2)
  )
  



# Prelim analyses ---------------------------------------------------------


dat_mds_reg <- dat_mds %>% 
  
  # Exclude 86 patients with missing relapse info (we will revisit this)
  filter(!is.na(ci_allo1)) %>% 
  
  # Admin censoring, time in years
  mutate(
    ci_allo1 = ci_allo1 / 12,
    ci_allo1 = ifelse(ci_allo1 >= 10, 10, ci_allo1),
    ci_s_allo1 = ifelse(ci_allo1 == 10, 0, as.numeric(ci_s_allo1) - 1)
  ) %>% 
  
  # Group categories of t-cell depletion
  mutate(
    tceldepl_allo1 = case_when(
      tceldepl_allo1 == "yes invivo, no exvivo" ~ "invivo_only",
      tceldepl_allo1 == "yes invivo+exvivo" ~ "exvivo",
      tceldepl_allo1 == "yes exvivo, no invivo" ~ "exvivo",
      is.na(tceldepl_allo1) ~ NA_character_,
      TRUE ~ "no" 
    ),
    tceldepl_allo1 = factor(
      tceldepl_allo1, 
      levels = c("no", "invivo_only", "exvivo")
    )
  ) %>% 
  
  # Make Karnofsky into ordered three cats
  mutate(
    KARNOFSK_allo1 = cut(
      x = KARNOFSK_allo1, 
      breaks = c(-Inf, 70, 80, 100), 
      include.lowest = T, ordered_result = T, 
      labels = c("<=70", "80", ">=90")
    )
  ) %>% 
  
  # Drop levels of factors class
  mutate_if(is.factor, ~ droplevels(.)) %>% 
  
  # Re-order crnocr
  mutate(
    crnocr = factor(crnocr, 
                    levels = c("Untreated/not aimed at remission",
                               "noCR", "CR"))
  )


# Table 1 -----------------------------------------------------------------

library(gtsummary)

dat_mds_reg %>% 
  select(-CENTRE_allo1, -srv_s_allo1, -srv_allo1,
         -AA_AUTO, -pb_diag1, -pb_allo1, -bm_diag1,
         -ci_allo1) %>% 
  tbl_summary(by = mdsclass) %>%
  add_n(statistic = "{p_miss}%", col_label = "% missing")

# Specifically for missingness
naniar::gg_miss_var(dat_mds_reg, facet = mdsclass, show_pct = T)

# Density plots MDS / cross tables with CR --------------------------------


blasts <- stringr::str_subset(names(dat_mds_reg), "bm_|pb_")


list_densplots <- purrr::map(
  .x = blasts,
  .f = ~ {
    dat_mds %>% 
      ggplot(aes(x = .[, .x], fill = mdsclass)) +
      geom_density(alpha= .5, adjust = 1) +
      xlim(c(0, 50)) +
      ggtitle(.x) +
      theme(axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) 
  }
)


library(ggpubr)

p <- ggarrange(plotlist = list_densplots, nrow = 2, ncol = 2, 
               common.legend = T, legend = "top")

annotate_figure(p, left = "Density", bottom = "Blast percentage")


# For CR groups
blasts <- stringr::str_subset(names(dat_mds_reg), "bm_allo1|pb_allo1")


list_densplots <- purrr::map(
  .x = blasts,
  .f = ~ {
    dat_mds %>% 
      ggplot(aes(x = .[, .x], fill = crnocr)) +
      geom_density(alpha= .5, adjust = 1, na.rm = T) +
      xlim(c(0, 50)) +
      ggtitle(.x) +
      theme(axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) 
  }
)


library(ggpubr)

p <- ggarrange(plotlist = list_densplots,# nrow = 2, 
               ncol = 2, 
               common.legend = T, legend = "top")

annotate_figure(p, left = "Density", bottom = "Blast percentage")


# Take maximum
dat_mds_reg %>% 
  mutate(max_blast = pmax(pb_allo1, pb_diag1, na.rm = T)) %>% 
  select(max_blast, mdsclass) %>% 
  tbl_summary(by = mdsclass, st)


# Cross tables crnocr and mds class
dat_mds_reg %>% 
  select(crnocr, mdsclass) %>% 
  tbl_summary(by = mdsclass)

dat_mds_reg %>% 
  select(crnocr, bm_allo1) %>% 
  tbl_summary(by = crnocr)


# Univ cox models relapse -------------------------------------------------



#  Vector of possible predictors 
predictors <- c("match_allo1_1", "age_allo1", "YEAR_allo1",
                "mdsclass", "donorrel", "bm_allo1", "AGEDONOR_allo1_1",
                "KARNOFSK_allo1", "crnocr", "tceldepl_allo1")
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
summary(univ_mods_rel$KARNOFSK_allo1)
plot(cox.zph(univ_mods_rel$KARNOFSK_allo1))

# Or all
purrr::map(univ_mods_rel, summary)

# Plot all coxzphs
purrr::map(univ_mods_rel, ~ plot(cox.zph(.x)))



# Univ cox models nrm -----------------------------------------------------


# Univariate cause-specific models (relapse)
univ_mods_nrm <- purrr::map(
  .x = predictors,
  .f = ~ {
    form <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", .x))
    mod <- coxph(form, data = dat_mds_reg)
  }
)

# Print any summary
summary(univ_mods_nrm$KARNOFSK_allo1)
plot(cox.zph(univ_mods_nrm$KARNOFSK_allo1))

# Or all
purrr::map(univ_mods_nrm, summary)

# Plot all coxzphs
purrr::map(univ_mods_nrm, ~ plot(cox.zph(.x)))




# Multivariable cox models ------------------------------------------------


preds <- paste(predictors, collapse = " + ")
preds 
  

# For relapse
form_rel <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 1) ~ ", preds))
form_rel
mod_rel <- coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ PATSEX + age_allo1 + 
                   YEAR_allo1 + mdsclass + donorrel + bm_allo1 + AGEDONOR_allo1_1 + 
                   KARNOFSK_allo1 + crnocr, 
                 data = dat_mds_reg)

mod_rel <- coxph(formula = Surv(ci_allo1, ci_s_allo1 == 1) ~ PATSEX + age_allo1 + 
                   YEAR_allo1 + mdsclass + donorrel + #bm_allo1 + 
                   AGEDONOR_allo1_1 + 
                   KARNOFSK_allo1 + crnocr, 
                 data = dat_mds_reg)

summary(mod_rel)
cox.zph(mod_rel)

cormat_rel <- cov2cor(mod_rel$var)
colnames(cormat_rel) <- rownames(cormat_rel) <- names(mod_rel$coefficients)
cormat_rel %>%View()


# For nrm
form_nrm <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", preds))
form_nrm
mod_nrm <- coxph(formula = Surv(ci_allo1, ci_s_allo1 == 2) ~ PATSEX + age_allo1 + 
                   YEAR_allo1 + mdsclass + donorrel + #bm_allo1 + 
                   AGEDONOR_allo1_1 + 
                   KARNOFSK_allo1 + crnocr, data = dat_mds_reg)

summary(mod_nrm)
cox.zph(mod_nrm)





# Prep mice matrix --------------------------------------------------------


# Which vars actually have some data missing?
var_names_miss <- naniar::miss_var_which(dat_prepped)

# Set methods accordingly
meth_miss <- setNames(
  c("logreg", rep("norm", 5), "", "logreg", "logreg", "polyreg"), 
 var_names_miss
)

meth_miss

# Prep imputation methods
meths <- setNames(rep("", ncol(dat_prepped)), names(dat_prepped))
meths[var_names_miss] <- meth_miss

meths

# Create matrix for MI
matpred <- matrix(1, ncol(dat_prepped), ncol(dat_prepped),
                  dimnames = list(names(dat_prepped), names(dat_prepped)))

diag(matpred) <- 0 # dont impute a var using itself

# We don't use time in imputation models
matpred[, c("time", "delta", "CENTRE_allo1", "rel_s_allo1",
            "srv_allo1", "srv_s_allo1", "ci_allo1", "ci_s_allo1",
            "rel_allo1", "AA_AUTO")] <- 0

matpred

# We only impute donor age
matpred[!(rownames(matpred) %in% var_names_miss), ] <- 0
matpred["rel_allo1",] <- 0




# Run mice ----------------------------------------------------------------


imps_mice <- mice(
  data = dat_prepped, 
  m = 5, 
  method = meths, 
  predictorMatrix = matpred
)



# Run smcfcs --------------------------------------------------------------

meths["cmv_combi_allo1_1"] <- "mlogit"

# you need the == ==
# Remove rel 
imps_smcfcs <- smcfcs::smcfcs(
  originaldata = dat_prepped, 
  smtype = "compet", 
  smformula = c(
    "Surv(time, delta == 1) ~ age_allo1 + AGEDONOR_allo1_1 + 
                bm_allo1 + ric_allo1 + mdsclass + donorrel +
                PATSEX + bm_allo1 + pb_allo1 + karnofsky + tbi_allo1 + cmv_combi_allo1_1",
    "Surv(time, delta == 2) ~ age_allo1 + AGEDONOR_allo1_1 + 
                bm_allo1 + ric_allo1 + mdsclass + donorrel +
                PATSEX + bm_allo1 + pb_allo1 + karnofsky + tbi_allo1 + cmv_combi_allo1_1"
  ),
  m = 5,
  method = meths
)





imps_smcfcs <- smcfcs::smcfcs(
  originaldata = dat_prepped, 
  smtype = "compet", 
  smformula = c(
    "Surv(time, delta == 1) ~ 1",
    "Surv(time, delta == 2) ~ 1"
  ),
  m = 5,
  method = meths
)
