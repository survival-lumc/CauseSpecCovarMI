
dat_mds <- fst::read_fst("data/dat-mds_admin-cens.fst") %>% data.table::setDT()



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


# Label factor levels -----------------------------------------------------


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


# Merge -------------------------------------------------------------------

tot <- merge(levels_dat, vars_meta, by = "var_name", all.y = TRUE)

# Add counts per factor
counters <- lapply(names(dat_mds), function(col) {
  
  if (is.numeric(dat_mds[[col]])) {
    counts <- dat_mds[, .(
      levels = NA_character_,
      count = sum(!is.na(get(col)))
    )]
    return(counts)
  } else {
    counts <- na.omit(dat_mds[, .(count = .N), by = col])
    data.table::setnames(counts, col, "levels")
    
    return(counts[, levels := as.character(levels)])
  }
})

names(counters) <- names(dat_mds)

final <- data.table::rbindlist(counters, idcol = "var_name") %>% 
  merge(tot, by = c("levels", "var_name"), all.y = TRUE)

saveRDS(final, file = "data-dictionary.rds")


# Try table ---------------------------------------------------------------

final <- readRDS("data/data-dictionary.rds")
data.table::setorder(final, "var_name", "level_num")

library(kableExtra)
final[is.na(levels_lab), levels_lab := var_label]
final[, .N, keyby = var_label]

options(knitr.table.format = "latex")

data.table::setnames(
  final,
  c("levels_lab", "count", "var_description"),
  c("Variable", "Count", "Description")
)



# One option
kbl(
  x = final[!(var_name  %in% c("srv_s_allo1", "srv_allo1")), 
            c("Variable", "Count", "Description")],
  #format = "latex",
  booktabs = "T", 
  linesep = ""
) %>% 
  kable_styling(full_width = T) %>% 
  column_spec(3, width = "22em") %>% 
  column_spec(2, width = "5em") %>% 
  collapse_rows(3, latex_hline = "none", valign = "top") %>% 
  # More programmatically here
  pack_rows(
    bold = F,
    index = c(
      " " = 4,
      "CMV Pat/Don" = 4,
      "CR stage" = 3,
      "Cytogenetics" = 3,
      "Don relation" = 2,
      "Comorbidity score" = 3,
      "Karnofsky" = 3,
      "Pat/Don sex match" = 4,
      "MDS class" = 3,
      " " = 2
    )
  ) 
   #%>%
  #column_spec(2, width = "5em")



# Other option, no counts

data.table::setnames(
  final,
  c("Variable", "var_label", "prop_miss"),
  c("Levels", "Variable", "% Missing")
)

final[Variable == Levels, Levels := ""]


kbl(
  x = final[!(var_name  %in% c("srv_s_allo1", "srv_allo1")), 
            c("Variable", "Description", "Levels", "% Missing")],
  format = "latex",
  booktabs = "T", 
  linesep = ""#,
  #escape = F
) %>% 
  kable_styling(font_size = 7) %>% #full_width = T) %>% 
  #kable_styling(latex_options = "scale_down") %>% 
  column_spec(2, width = "20em") %>% 
  collapse_rows(2, latex_hline = "none", valign = "top") %>%
  collapse_rows(1, latex_hline = "none", valign = "top") %>% 
  collapse_rows(4, latex_hline = "none", valign = "top")
  
# Maybe with type columns?

