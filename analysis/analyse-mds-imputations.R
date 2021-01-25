##***************************************##
## Compare synthetic and mds imputations ##
##***************************************##


devtools::load_all()

# Read-in
dat_mds <- fst::read_fst("data/dat-mds_admin-cens.fst") %>% 
  data.table::setDT()

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2)) 


# Prep formulas -----------------------------------------------------------



outcomes <- c("ci_s_allo1", "ci_allo1", "srv_s_allo1", "srv_allo1")
predictors <- sort(colnames(dat_mds)[!(colnames(dat_mds) %in% outcomes)]) 

# Both REL and NRM have same rhs
rhs <- paste(predictors, collapse = " + ")

# Make both model formulas
form_rel <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 1) ~ ", rhs))
form_nrm <- as.formula(paste0("Surv(ci_allo1, ci_s_allo1 == 2) ~ ", rhs))



# Load imputation ---------------------------------------------------------

# Load in objects
imps_mice <- readRDS("data/imps-mds-mice.rds")
imps_mice_synth <- readRDS("data/imps-mds-mice-synth.rds")
imps_smcfcs <- readRDS("data/imps-mds-smcfcs.rds")
imps_smcfcs_synth <- readRDS("data/imps-mds-smcfcs-synth.rds")

# Plot the synth ones
plot(imps_mice_synth)
plot(imps_smcfcs_synth)

# Prepare lists 
impdats_mice <- mice::complete(imps_mice, action = "all")
impdats_smcfcs <- imps_smcfcs$impDatasets
impdats_mice_synth <- mice::complete(imps_mice_synth, action = "all")
impdats_smcfcs_synth <- imps_smcfcs_synth$impDatasets


# Relapse models ----------------------------------------------------------

mice_rel <- lapply(
  impdats_mice, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)

# Mice models
mice_rel_synth <- lapply(
  impdats_mice_synth, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)

# smcfcs
smcfcs_rel <- lapply(
  impdats_smcfcs, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)

smcfcs_rel_synth <- lapply(
  impdats_smcfcs_synth, 
  function(imp) survival::coxph(form_rel, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)


# NRM models --------------------------------------------------------------


mice_nrm <- lapply(
  impdats_mice, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)

# Mice models
mice_nrm_synth <- lapply(
  impdats_mice_synth, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)

# smcfcs
smcfcs_nrm <- lapply(
  impdats_smcfcs, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)

smcfcs_nrm_synth <- lapply(
  impdats_smcfcs_synth, 
  function(imp) survival::coxph(form_nrm, data = imp)
) %>% 
  mice::pool() %>% 
  summary(conf.int = TRUE)


# Comparison --------------------------------------------------------------


# Mice relapse
cbind.data.frame(
  "var_rel" = mice_rel$term, 
  "orig" = mice_rel$estimate, 
  "synth" = mice_rel_synth$estimate
)

# Mice nrm
cbind.data.frame(
  "var_nrm" = mice_nrm$term, 
  "orig" = mice_nrm$estimate, 
  "synth" = mice_nrm_synth$estimate
)


# smcfcs relapse
cbind.data.frame(
  "var_rel" = smcfcs_rel$term, 
  "orig" = smcfcs_rel$estimate, 
  "synth" = smcfcs_rel_synth$estimate
)

# smcfcs nrm
cbind.data.frame(
  "var_nrm" = smcfcs_nrm$term, 
  "orig" = smcfcs_nrm$estimate, 
  "synth" = smcfcs_nrm_synth$estimate
)



# CCA ---------------------------------------------------------------------


dat_mds_synth <- fst::read_fst("data/dat-mds-synth.fst") %>% 
  data.table::setDT()


mod_rel <- survival::coxph(form_rel, data = dat_mds)
mod_rel_synth <- survival::coxph(form_rel, data = dat_mds_synth)
mod_nrm <- survival::coxph(form_nrm, data = dat_mds)
mod_nrm_synth <- survival::coxph(form_nrm, data = dat_mds_synth)


# Relapse
cbind.data.frame(
  "orig" = coef(mod_rel),
  "synth" = coef(mod_rel_synth)
)

# NRM
cbind.data.frame(
  "orig" = coef(mod_nrm),
  "synth" = coef(mod_nrm_synth)
)


