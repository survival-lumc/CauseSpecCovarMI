##************************##
## Process/clean MDS data ##
##************************##


# Reading-in --------------------------------------------------------------


# Read-in data with hctci
dat_hctci <- sjlabelled::read_spss(
  path = "data-raw/Final_MDS_Longterm_20200903_LK.sav" 
) %>% 
  
  # Clean variable names (i.e. bring lower case, with _ in between)
  janitor::clean_names() %>% 
  
  # Make into data.table
  data.table::setDT() %>% 
  
  # Subset aa_auto for matching and new vars
  .[, .(
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
dat_mds <- data.table::copy(dat_combined) %>% 
  
  # Select variables, Remove one patient with wrong donor age (data error)
  .[aa_auto != 357439, ..vars_keep] %>% 
  
  # Set time in years, age and donor age in decades
  .[, ':=' (
    ci_allo1 = ci_allo1 / 12,
    srv_allo1 = srv_allo1 / 12,
    agedonor_allo1_decades = agedonor_allo1_1 / 10,
    age_allo1_decades = age_allo1 / 10
  )] %>% 
  
  # For some vars, use the spss values labels as factor labels
  .[, (vars_val_labs) := lapply(.SD, sjlabelled::as_label), 
    .SDcols = vars_val_labs] %>% 
  
  # Make Karnofsky into ordered three categories
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
  
  # Reorder crnocr, order hctci
  .[, ':=' (
    crnocr = factor(
      crnocr, 
      levels = c("CR", "noCR", "Untreated/not aimed at remission")
    ),
    hctci_risk = as.ordered(hctci_risk)
  )] %>%
  
  # Delete redundant age and donor age not in decades (they are in dat_orig)
  .[, ':=' (
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
fst::write_fst(x = dat_mds, path = "analysis/data/dat-mds.fst")

# Save a second version with admin censoring at 10 years
dat_mds[, ':=' (
  ci_allo1 = ifelse(ci_allo1 >= 10, 10, ci_allo1),
  srv_allo1 = ifelse(srv_allo1 >= 10, 10, srv_allo1)
)] 

dat_mds[, ':=' (
  ci_s_allo1 = ifelse(ci_allo1 == 10, 0, sjlabelled::as_numeric(ci_s_allo1)),
  srv_s_allo1 = ifelse(srv_allo1 == 10, 0, sjlabelled::as_numeric(srv_s_allo1))
)]

fst::write_fst(x = dat_mds, path = "analysis/data/dat-mds_admin-cens.fst")



# Create data dictionary --------------------------------------------------


# Keep original names
vars <- setNames(rep(list(""), ncol(dat_mds)), names(dat_mds))

# Add the descriptions
vars[["ci_s_allo1"]] <- cbind.data.frame(
  "var_label" = "Competing event indicator",
  "var_description" = "Competing event indicator"
)

vars[["ci_allo1"]] <- cbind.data.frame(
  "var_label" = "Time to competing event",
  "var_description" = "Time from alloHCT to competing event (months)"
)

vars[["srv_s_allo1"]] <- cbind.data.frame(
  "var_label" = "Death/censoring indicator",
  "var_description" = "Time from alloHCT to competing event (years)"
)

vars[["srv_allo1"]] <- cbind.data.frame(
  "var_label" = "Time to death/censoring",
  "var_description" = "Time from alloHCT to death or censoring (years)"
)

vars[["srv_allo1"]] <- cbind.data.frame(
  "var_label" = "Time to death/censoring",
  "var_description" = "Time from alloHCT to death or censoring (years)"
)

vars[["match_allo1_1"]] <- cbind.data.frame(
  "var_label" = "Pat/Don sex match",
  "var_description" = "Sex match patient and donor"
)

vars[["mdsclass"]] <- cbind.data.frame(
  "var_label" = "MDS class",
  "var_description" = "MDS groups based on subclassification (WHO and FAB) at alloHCT"
)

vars[["donorrel"]] <- cbind.data.frame(
  "var_label" = "Don relation",
  "var_description" = "Relation with donor"
)

vars[["karnofsk_allo1"]] <- cbind.data.frame(
  "var_label" = "Karnofsky",
  "var_description" = "Karnofsky or Lansky status"
)

vars[["crnocr"]] <- cbind.data.frame(
  "var_label" = "CR stage",
  "var_description" = "Stage at alloHCT"
)

vars[["cmv_combi_allo1_1"]] <- cbind.data.frame(
  "var_label" = "CMV Pat/Don",
  "var_description" = "Cytomegalovirus in patient and donor"
)

vars[["cytog_threecat"]] <- cbind.data.frame(
  "var_label" = "Cytogenetics",
  "var_description" = "Cytogenetics score"
)

vars[["hctci_risk"]] <- cbind.data.frame(
  "var_label" = "Comorbidity score",
  "var_description" = "HCTCI comorbidity index score"
)

vars[["agedonor_allo1_decades"]] <- cbind.data.frame(
  "var_label" = "Age (Don)",
  "var_description" = "Donor age at alloHCT (decades)"
)

vars[["age_allo1_decades"]] <- cbind.data.frame(
  "var_label" = "Age (Pat)",
  "var_description" = "Patient age at alloHCT (decades)"
)

missing_dat <- data.table::transpose(
  dat_mds[, lapply(.SD, function(col) mean(is.na(col)))], 
  keep.names = "var_name"
)
data.table::setnames(missing_dat, setdiff(names(missing_dat), "var_name"), "prop_miss")
vars_meta <- merge(data.table::rbindlist(vars, idcol = "var_name"), missing_dat)
vars_meta[, prop_miss := paste0(round(prop_miss, 4) * 100)]


# Do the same for factor levels 
factors <- sapply(dat_mds, is.factor)
levs <- lapply(names(dat_mds)[factors], function(col) {
  var <- dat_mds[[col]]
  cbind.data.frame("levels" = levels(var), "level_num" = 1:length(levels(var)))
})
names(levs) <- names(dat_mds)[factors]

levs[["match_allo1_1"]] <- transform(
  levs[["match_allo1_1"]],
  "levels_lab" = c("M/M", "M/F", "F/M", "F/F")
)

levs[["mdsclass"]] <- transform(
  levs[["mdsclass"]],
  "levels_lab" = c("MDS w/o excess blasts", "MDS w/ excess blasts", "sAML")
)

levs[["donorrel"]] <- transform(
  levs[["donorrel"]],
  "levels_lab" = c("Identical sibling", "Other")
)

levs[["karnofsk_allo1"]] <- transform(
  levs[["karnofsk_allo1"]],
  "levels_lab" = c(">=90", "80", "<=70")
)

levs[["crnocr"]] <- transform(
  levs[["crnocr"]],
  "levels_lab" = c("CR", "no CR", "Untreated")
)

levs[["cmv_combi_allo1_1"]] <- transform(
  levs[["cmv_combi_allo1_1"]],
  "levels_lab" = c("-/-", "-/+", "+/-", "+/+")
)

levs[["cytog_threecat"]] <- transform(
  levs[["cytog_threecat"]],
  "levels_lab" = c("Good (<=3)", "Poor (4)", "Very poor (5)")
)

levs[["hctci_risk"]] <- transform(
  levs[["hctci_risk"]],
  "levels_lab" = c("Low risk (0)", "Interm. risk (1-2)", "High risk (>=3)")
)

levels_dat <- data.table::rbindlist(levs, idcol = "var_name")


vlapply(names(dato), function(col) {
  
  if (is.numeric(dato[[col]])) {
    
    summ <- data.table::data.table(
      level = NA,
      count
    )
    #summ[, pct := count / sum(count)]
    #data.table::setnames(summ, old = nam, "levs")
  } else {
    
    summ <- dato[, .(count = .N), keyby = col]
    summ[, pct := round(count / sum(count), 2)]
    levs <- levels(dato[[col]])
    data.table::set(summ, j = "level_num", value = match(summ[[col]], levs))
    data.table::setnames(summ, col, "level")
  }
  
  return(summ)
})

