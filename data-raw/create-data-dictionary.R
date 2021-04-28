##************************##
## Create data dictionary ##
##************************##

# Choose whether synthetic or not
synth <- FALSE

if (synth) {
  dat_mds <- CauseSpecCovarMI::dat_mds_synth %>% data.table::setDT()
} else {
  dat_mds <- fst::read_fst("data/dat-mds_admin-cens.fst") %>% data.table::setDT()
  dat_mds[, c("srv_s_allo1", "srv_allo1") := NULL]
}


# Variable descriptions ---------------------------------------------------


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

vars[["match_allo1_1"]] <- cbind.data.frame(
  "var_label" = "Patient/Donor sex match",
  "var_description" = "Sex match patient and donor"
)

vars[["mdsclass"]] <- cbind.data.frame(
  "var_label" = "MDS class",
  "var_description" = "MDS groups based on subclassification at alloHCT"
)

vars[["donorrel"]] <- cbind.data.frame(
  "var_label" = "HLA match patient/donor",
  "var_description" = "HLA match between patient and donor"
)

vars[["karnofsk_allo1"]] <- cbind.data.frame(
  "var_label" = "Karnofsky",
  "var_description" = "Karnofsky performance status"
)

vars[["crnocr"]] <- cbind.data.frame(
  "var_label" = "Stage",
  "var_description" = "Stage at alloHCT"
)

vars[["cmv_combi_allo1_1"]] <- cbind.data.frame(
  "var_label" = "CMV Patient/Donor",
  "var_description" = "CMV status in patient and donor"
)

vars[["cytog_threecat"]] <- cbind.data.frame(
  "var_label" = "Cytogenetics",
  "var_description" = "Cytogenetics categories used for IPSS-R"
)

vars[["hctci_risk"]] <- cbind.data.frame(
  "var_label" = "Comorbidity score",
  "var_description" = "HCT-CI score"
)

vars[["agedonor_allo1_decades"]] <- cbind.data.frame(
  "var_label" = "Age (Donor)",
  "var_description" = "Donor age at alloHCT (decades)"
)

vars[["age_allo1_decades"]] <- cbind.data.frame(
  "var_label" = "Age (Patient)",
  "var_description" = "Patient age at alloHCT (decades)"
)

# Keep proportion of missing data per variable
missing_dat <- data.table::transpose(
  dat_mds[, lapply(.SD, function(col) mean(is.na(col)))], 
  keep.names = "var_name"
)

data.table::setnames(missing_dat, setdiff(names(missing_dat), "var_name"), "prop_miss")
vars_meta <- merge(data.table::rbindlist(vars, idcol = "var_name"), missing_dat)
vars_meta[, prop_miss := paste0(round(prop_miss, 4) * 100)]


# Label factor levels -----------------------------------------------------


# Extract factors and their levels in a list
factors <- sapply(dat_mds, is.factor)
levs <- lapply(names(dat_mds)[factors], function(col) {
  var <- dat_mds[[col]]
  cbind.data.frame("levels" = levels(var), "level_num" = seq_len(length(levels(var))))
})

names(levs) <- names(dat_mds)[factors]


# Add factor level labels (for data dictionary table)
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
  "levels_lab" = c("HLA-identical sibling", "Other")
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
  "levels_lab" = c("V. good/good/interm.", "Poor", "V. poor")
)

levs[["hctci_risk"]] <- transform(
  levs[["hctci_risk"]],
  "levels_lab" = c("Low risk (0)", "Interm. risk (1-2)", "High risk (>=3)")
)

levels_dat <- data.table::rbindlist(levs, idcol = "var_name")


# Keep all label info for vars and factor in one df -----------------------


labels_all <- merge(levels_dat, vars_meta, by = "var_name", all.y = TRUE)

# Add counts per factor - and REL / NRM counts too
counters <- lapply(names(dat_mds), function(col) {
  
  if (is.numeric(dat_mds[[col]])) {
    counts <- dat_mds[, .(
      levels = NA_character_,
      count = .N,
      count_REL = sum(ci_s_allo1 == 1),
      count_NRM = sum(ci_s_allo1 == 2)
    ), by = is.na(get(col))]
    data.table::setnames(counts, "is.na", "miss_ind")
    return(counts[miss_ind == FALSE])
  } else {
    counts <- na.omit(dat_mds[, .(
      count = .N,
      count_REL = sum(ci_s_allo1 == 1),
      count_NRM = sum(ci_s_allo1 == 2)
    ), by = col])
    data.table::setnames(counts, col, "levels")
    
    return(counts[, levels := as.character(levels)])
  }
})

names(counters) <- names(dat_mds)

# Bind counters, labels all together
dictionary_df <- data.table::rbindlist(
  l = counters, 
  idcol = "var_name", 
  fill = TRUE
) %>% 
  data.table::merge.data.table(
    y = labels_all, 
    by = c("levels", "var_name"), 
    all.y = TRUE
  )

dictionary_df[, "miss_ind" := NULL]

# Save if not the synthetic version - as "internal"data
if (!synth) save(dictionary_df, file = "R/data_dictionary.rda")

# Descriptives table ------------------------------------------------------


dictionary_df <- get(load("R/data_dictionary.rda"))
data.table::setorder(dictionary_df, "var_name", "level_num")
dictionary_df[is.na(levels_lab), levels_lab := ""]

# Make columns for LaTex
dictionary_df[, levels_lab_tex := data.table::fcase(
  var_name == "cmv_combi_allo1_1", paste0("$", levels_lab, "$"),
  grepl(pattern = ">=", x = levels_lab), gsub(">=", x = levels_lab, replacement = "$\\\\geq$"),
  grepl(pattern = "<=", x = levels_lab), gsub("<=", x = levels_lab, replacement = "$\\\\leq$")
)]

dictionary_df[is.na(levels_lab_tex), levels_lab_tex := levels_lab]

# Set column names
data.table::setnames(
  dictionary_df,
  c("var_label", "var_description", "levels_lab_tex", "prop_miss"),
  c("Variable", "Description", "Levels", "\\% Missing")
)


caption <- "Data dictionary with predictor variables and their descriptions, levels and proportion missing data. Abbrevations: CMV = cytomegalovirus, CR = complete remission, IPSS-R = International Prognostic Scoring System, V. = very, interm. = intermediate, HLA = Human leukocyte antigen, HCT-CI = Hematopoietic stemcell transplantation-comorbidity index, M = male, F = female, MDS = myelodysplastic syndromes, sAML = secondary acute myeloid leukemia, w/ = with, w/o = without."

# Make table
dictionary_df[!(var_name  %in% c("srv_s_allo1", "srv_allo1", "ci_allo1", "ci_s_allo1")), c(
  "Variable", "Description", "Levels", "\\% Missing"
)] %>% 
  kableExtra::kbl(
    format = "latex",
    booktabs = "T", 
    position = "h",
    caption = caption,
    linesep = "",
    escape = F, 
    digits = 2
  ) %>% 
  kableExtra::kable_styling(font_size = 7) %>% #full_width = T) %>% 
  kableExtra::column_spec(2, width = "20em") %>% 
  kableExtra::collapse_rows(2, latex_hline = "none", valign = "top") %>%
  kableExtra::collapse_rows(1, latex_hline = "none", valign = "top") %>% 
  kableExtra::collapse_rows(4, latex_hline = "none", valign = "top")


# Table with extra measurment level column --------------------------------


classes <- vapply(dat_mds, function(col) {
  cls <- class(col)[1]
  if (length(levels(col)) == 2) cls <- "Binary"
  return(cls)
}, FUN.VALUE = character(1))

dictionary_df[, "Meas. level" := classes[match(
  dictionary_df$var_name, 
  names(classes)
)]]

dictionary_df[, "Meas. level" := data.table::fcase(
  `Meas. level` == "numeric", "Continuous",
  `Meas. level` == "factor", "Nominal categorical",
  `Meas. level` == "ordered", "Ordered categorical",
  `Meas. level` == "Binary", "Binary"
)]

# Edit for mds class, which is ordered
dictionary_df[var_name == "mdsclass", `Meas. level` := "Ordered categorical"]

dictionary_df[!(is.na(level_num) | level_num == 1), `Meas. level` := ""]


dictionary_df <- dictionary_df[!(var_name  %in% c("srv_s_allo1", "srv_allo1", "ci_allo1", "ci_s_allo1")), c(
  "Variable", "Description", "Meas. level", "Levels", "\\% Missing"
)]

dictionary_df %>% 
  kableExtra::kbl(
    format = "latex",
    booktabs = "T", 
    position = "h",
    caption = caption,
    linesep = "",
    escape = F, 
    digits = 2
  ) %>% 
  kableExtra::kable_styling(font_size = 8) %>% 
  kableExtra::column_spec(1, width = "7em") %>% 
  kableExtra::column_spec(2, width = "10em") %>% 
  kableExtra::collapse_rows(1, latex_hline = "none", valign = "top") %>% 
  kableExtra::collapse_rows(2, latex_hline = "none", valign = "top") %>% 
  kableExtra::collapse_rows(5, latex_hline = "none", valign = "top")


