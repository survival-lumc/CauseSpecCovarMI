

# After reading-in imps MICE
imp_dats <- mice::complete(imps_mice, action = "all")

dato <- imp_dats[[1]]




tmat <- mstate::trans.comprisk(2, c("Rel", "NRM"))
covs <- predictors_mod1

# Long format
dat_msprepped <- mstate::msprep(time = c(NA, "ci_allo1", "ci_allo1"),
                                status = c(NA, "ev1", "ev2"), 
                                data = dato,
                                trans = tmat,
                                keep = covs) 

# Expand covariates - deal with variable names later with longnames
dat_expanded <- mstate::expand.covs(
  dat_msprepped, covs, append = TRUE, longnames = F
)

# Full mod
namos <- colnames(dat_expanded)
gsub(x = namos, pattern = "*(.[1-9+])", replacement = "-HELLO")
# use grep
#cond <-
rhs_mstate <- paste(
  names(data.table(dat_expanded)[, !(id:status)]),
  collapse = " + "
)

form_mstate <- as.formula(paste0("Surv(time, status) ~ ", rhs_mstate))


# Function generating reference patient 
survival::coxph(form_mstate, dat_expanded)


#  Attempt with riskregression --------------------------------------------

# USe median
lapply(dat_mds_reg, function(x) if (is.factor(x)) table(x) else median(x, na.rm = T))

imps_mice[[1]]

library(riskRegression)

riskreg_form <- as.formula(paste0("Hist(ci_allo1, ci_s_allo1) ~ ", preds))

mod_test <- CSC(riskreg_form, data = dato)

ref_pat <- dat_mds_reg[1, 7:ncol(dat_mds_reg)]

newvals <- c("recipient male - donor male", "sAML", "Identical sibling",
                  ">=90", "CR", "-/-", "good(<=3)", "low risk (0)",
                  4.245213, 5.331636)
nmz <- colnames(ref_pat)
ref_pat[1, names(ref_pat) := as.list(newvals)]
ref_pat

refpato <- rbind(ref_pat, ref_pat, ref_pat)
refpato[, mdsclass := c("MDS without excess blasts", 
                        "MDS with excess blasts",
                        "sAML")]


ev1 <- print(predict(mod_test, newdata = refpato, cause = 1, time = 5, se = T, confint = F))
ev2 <- print(predict(mod_test, newdata = refpato, cause = 2, time = 5, se = T, confint = F))

res <- rbind(ev1, ev2, idcol = "cause")[, c("cause", "mdsclass", "absRisk", "absRisk.se")]

pred_impos <- function(impdat) {
  
  modo <- CSC(riskreg_form, data = impdat)
  ev1 <- print(predict(modo, newdata = refpato, cause = 1, 
                       time = 5, se = T, confint = F))
  ev2 <- print(predict(modo, newdata = refpato, cause = 2,
                       time = 5, se = T, confint = F))
  
  res <- rbind(ev1, ev2, idcol = "cause")[, c("cause", "mdsclass", 
                                              "absRisk", "absRisk.se")]
  
  return(res)
}

cloglog <- Vectorize(function(x) log(-log(1 - x)))
inv_cloglog <- Vectorize(function(x) 1 - exp(-exp(x)))


mice_res_preds <- lapply(imp_dats, pred_impos)


testo <- rbindlist(mice_res_preds) 

mice_predos <- testo[, ':=' (
    p  = absRisk,
    p_trans = cloglog(absRisk),
    var_p = (absRisk.se)^2
  )] %>% 
  .[, ':=' (
    Ui = var_p / (log(1 - p) * (1 - p))^2
  )] %>% 
  .[, .(
    m = .N,
    Ubar = mean(Ui),
    B = var(p_trans),
    Qbar = mean(p_trans)
  ), by = c("cause", "mdsclass")] %>% 
  .[, ':=' (
    total_var = Ubar + (1 + m^-1) * B
  ), by = c("cause", "mdsclass")] %>% 
  
  # Transf back
  .[, .(
    p_pooled = inv_cloglog(Qbar),
    #se_pooled = inv_cloglog(sqrt(total_var)),
    CI_low = inv_cloglog(Qbar - qnorm(0.975) * sqrt(total_var)),
    CI_upp = inv_cloglog(Qbar + qnorm(0.975) * sqrt(total_var))
  ), by = c("cause", "mdsclass")] %>% 
  print()



# smcfcs ------------------------------------------------------------------

implist_smcfcs


smcfcs_res_preds <- lapply(implist_smcfcs$imputations, pred_impos)


testo_smcfcs <- rbindlist(smcfcs_res_preds) 

smcfcs_predos <- testo_smcfcs[, ':=' (
  p  = absRisk,
  p_trans = cloglog(absRisk),
  var_p = (absRisk.se)^2
)] %>% 
  .[, ':=' (
    Ui = var_p / (log(1 - p) * (1 - p))^2
  )] %>% 
  .[, .(
    m = .N,
    Ubar = mean(Ui),
    B = var(p_trans),
    Qbar = mean(p_trans)
  ), by = c("cause", "mdsclass")] %>% 
  .[, ':=' (
    total_var = Ubar + (1 + m^-1) * B
  ), by = c("cause", "mdsclass")] %>% 
  
  # Transf back
  .[, .(
    p_pooled = inv_cloglog(Qbar),
    #se_pooled = inv_cloglog(sqrt(total_var)),
    CI_low = inv_cloglog(Qbar - qnorm(0.975) * sqrt(total_var)),
    CI_upp = inv_cloglog(Qbar + qnorm(0.975) * sqrt(total_var))
  ), by = c("cause", "mdsclass")] %>% 
  print()


# CCA ---------------------------------------------------------------------


modo_CCA <- CSC(riskreg_form, data = dat_mds_reg)
ev1 <- print(predict(modo_CCA, newdata = refpato, cause = 1, 
                     time = 5, se = T, confint = F))
ev2 <- print(predict(modo_CCA, newdata = refpato, cause = 2,
                     time = 5, se = T, confint = F))

res <- rbind(ev1, ev2, idcol = "cause")[, c("cause", "mdsclass", 
                                            "absRisk", "absRisk.se")]

# Also delta method for res
CCA_predos <- res %>% 
  .[, ':=' (
    p = absRisk,
    p_trans = cloglog(absRisk),
    var_p = absRisk.se^2
  )] %>% 
  .[, var_delta := var_p / (log(1 - p) * (1 - p))^2] %>% 
  .[, ':=' (
    CI_low = inv_cloglog(p_trans - qnorm(0.975) * sqrt(var_delta)),
    CI_upp = inv_cloglog(p_trans + qnorm(0.975) * sqrt(var_delta)),
    p_pooled = p
  )] %>% 
  .[, c("cause", "mdsclass", "p_pooled", "CI_low", "CI_upp")] %>% 
  print()


mice_predos
smcfcs_predos
res

# Do to 1 decimal!!
res_comb <- rbind(CCA_predos, mice_predos, smcfcs_predos, idcol = "method") %>% 
  .[, method := factor(method, levels = 1:3, labels = c("CC", "CH12", "smcfcs"))] %>% 
  .[, est := paste0(
    round(p_pooled * 100, 1),
    " [",
    round(CI_low * 100, 1),
    "; ",
    round(CI_upp * 100, 1),
    "]"
  )] %>% 
  dcast(formula = cause + mdsclass ~ method, value.var = "est") %>% 
  .[, mdsclass := factor(mdsclass, levels = c(
    "MDS without excess blasts",
    "MDS with excess blasts",
    "sAML"
  ))] %>% 
  .[order(cause, mdsclass)] %>% 
  .[, -1]


print(xtable(res_comb), include.rownames = F)
