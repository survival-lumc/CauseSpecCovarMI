##************************##
## Process/clean MDS data ##
##************************##

# Script to clean original MDS dataset - synthetic dataset already cleaned.

# Reading-in --------------------------------------------------------------


# Read-in data with hctci
dat_hctci <- sjlabelled::read_spss(
  path = "data-raw/Final_MDS_Longterm_20200903_LK.sav" 
) %>% 
  
  # Clean variable names (i.e. bring lower case, with _ in between)
  janitor::clean_names() %>% 
  data.table::setDT()
  
  # Subset aa_auto for matching and new vars
dat_hctci <- dat_hctci[, .(
  aa_auto,
  hctci,
  hctci_risk,
  agedonor_allo1_1
)]


# Load-in cytogenetics data
load("data-raw/cytodata.RData")
data.table::setDT(cytodata) 
data.table::setnames(cytodata, new = c("aa_auto", "cyto_score", "vchromos"))

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
  path = "data-raw/MDS_longterm_R_180614.sav"
) %>% 
  
  # Clean variable names
  janitor::clean_names() %>%
  data.table::setDT() 


# Combine all three files via Reduce + merge 
dat_combined <- Reduce(
  f = function(x, y) data.table::merge.data.table(x, y, all.x = T, by = "aa_auto"),
  x = list(dat_orig, dat_hctci, cytodata)
)

# Clean env
rm(dat_hctci, cytodata)


# Data / variable formatting ----------------------------------------------


# Define variables to keep - see table for motivation
vars_keep <- c(
  
  # Outcome vars
  "ci_s_allo1", 
  "ci_allo1", 
  "srv_s_allo1", 
  "srv_allo1",
  
  # aa_auto and (possible) covariates
  "aa_auto", 
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

# Factors that will use their SPSS value labels 
vars_val_labs <- c(
  "mdsclass",
  "crnocr",
  "donorrel",
  "hctci_risk",
  "cmv_combi_allo1_1",
  "match_allo1_1"
)

# Make dataset ready for regression
dat_mds <- data.table::copy(dat_combined)  
  
# Select variables, Remove one patient with wrong donor age (data error)
dat_mds <- dat_mds[aa_auto != 357439, ..vars_keep] 
  
# Set time in years, age and donor age in decades
dat_mds[, ':=' (
  ci_allo1 = ci_allo1 / 12,
  srv_allo1 = srv_allo1 / 12,
  agedonor_allo1_decades = agedonor_allo1_1 / 10,
  age_allo1_decades = age_allo1 / 10
)] 
  
# For some vars, use the spss values labels as factor labels
dat_mds[, (vars_val_labs) := lapply(.SD, sjlabelled::as_label), .SDcols = vars_val_labs] 
  
# Make Karnofsky into ordered three categories
dat_mds[, karnofsk_allo1 := cut(
  x = karnofsk_allo1, 
  breaks = c(-Inf, 70, 80, 100), 
  include.lowest = T, 
  ordered_result = T, 
  labels = c("<=70", "80", ">=90")
)] 
  
# Reverse ordering to have >= 90 as reference
dat_mds[, karnofsk_allo1 := factor(
  karnofsk_allo1, levels = rev(levels(karnofsk_allo1))
)]  
  
# Reorder crnocr, order hctci
dat_mds[, ':=' (
  crnocr = factor(
    crnocr, 
    levels = c("CR", "noCR", "Untreated/not aimed at remission")
  ),
  hctci_risk = as.ordered(hctci_risk)
)] 
  
# Delete redundant age and donor age not in decades (they are in dat_orig)
dat_mds[, ':=' (
  agedonor_allo1_1 = NULL,
  age_allo1 = NULL
)]  

# Drop unused factor levels 
which_factors <- names(dat_mds)[sapply(dat_mds, is.factor)]
dat_mds[, (which_factors) := lapply(.SD, droplevels), .SDcols = which_factors]

# Conserve var labels from original data (for mutated vars)
dat_mds <- data.table::copy(
  sjlabelled::copy_labels(df_new = dat_mds, df_origin = dat_orig)
)

# aa_auto no longer needed
dat_mds[, aa_auto := NULL]


# Export ------------------------------------------------------------------


# Save as fst
#fst::write_fst(x = dat_mds, path = "analysis/data/dat-mds.fst")

# Save a second version with admin censoring at 10 years
dat_mds[, ':=' (
  ci_allo1 = ifelse(ci_allo1 >= 10, 10, ci_allo1),
  srv_allo1 = ifelse(srv_allo1 >= 10, 10, srv_allo1)
)] 

dat_mds[, ':=' (
  ci_s_allo1 = ifelse(ci_allo1 == 10, 0, sjlabelled::as_numeric(ci_s_allo1)),
  srv_s_allo1 = ifelse(srv_allo1 == 10, 0, sjlabelled::as_numeric(srv_s_allo1))
)]

#fst::write_fst(x = dat_mds, path = "analysis/data/dat-mds_admin-cens.fst")

