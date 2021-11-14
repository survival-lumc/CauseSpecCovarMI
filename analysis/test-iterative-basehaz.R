devtools::load_all()
library(mice)
library(survival)

# https://github.com/lbeesleyBIOSTAT/SRMIMI_Example_Code/blob/main/Example_Code_Normal.R
# https://github.com/alexanderrobitzsch/miceadds
# https://www.gerkovink.com/miceVignettes/Passive_Post_processing/Passive_imputation_post_processing.html

baseline <- CauseSpecCovarMI::mds_shape_rates
shape_ev1 <- baseline[baseline$state == "REL", "shape"]
base_rate_ev1 <- baseline[baseline$state == "REL", "rate"]

ev1_pars <- list(
  "a1" = shape_ev1, 
  "h1_0" = base_rate_ev1,
  "b1" = .5, 
  "gamm1" = 1
)

# Parameters Weibull event 2
ev2_pars <- list(
  "a2" = baseline[baseline$state == "NRM", "shape"], 
  "h2_0" = baseline[baseline$state == "NRM", "rate"], 
  "b2" = .5, 
  "gamm2" = .5
)

dat <- generate_dat(
  n = 2000,
  X_type = "binary", 
  r = 0.5, 
  ev1_pars = ev1_pars,
  ev2_pars = ev2_pars, 
  rate_cens = baseline[baseline$state == "EFS", "rate"], 
  mech = "MAR", 
  p = 0.5,
  eta1 = -1
)


# Try here
mats <- get_predictor_mats(dat) 

# For binary mice = smcfcs
methods <- set_mi_methods(
  dat = dat, 
  var_names_miss = naniar::miss_var_which(dat), 
  imp_type = "mice", 
  cont_method = "norm" 
) 

m <- c(1) # Number of imputations of interest
iters_MI <- 1

# Try iterative updating with ch1: 
update_basehaz <- function(time, delta, x, z) {
  mod <- coxph(Surv(time, delta) ~ x + z, control = survival::coxph.control(timefix = FALSE))
  baseh_df <- basehaz(mod, centered = FALSE)
  haz <- baseh_df[match(time, baseh_df[["time"]]), ][["hazard"]]
  return(haz)
}

#update_basehaz(dat$t, dat$ev1, dat$X_orig, dat$Z)
methods["H1"] <- paste("~I(", expression(update_basehaz(t, ev1, X, Z)),")")
methods
predmat <- mats$CH1
predmat["H1", c("X", "Z", "t", "ev1", "ev2")] <- 1
predmat

dat$H1 <- NA
dat$H1

imp_ch1 <- mice::mice(
  dat, 
  m = m[length(m)],
  method = methods, 
  predictorMatrix = predmat,
  maxit = iters_MI#, 
  #print = FALSE
)

imp_ch1$imp$X

plot(imp_ch1)

#make.method

imps_comp <- mice::complete(imp_ch1, action = "all")
cbind.data.frame(imps_comp$`1`$H1,
                 imps_comp$`2`$H1,
                 imps_comp$`3`$H1) |>  View()
