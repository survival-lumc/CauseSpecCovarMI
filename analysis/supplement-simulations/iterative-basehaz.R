# References and load libraries -------------------------------------------

devtools::load_all()
library(mice)
library(survival)

# https://github.com/lbeesleyBIOSTAT/SRMIMI_Example_Code/blob/main/Example_Code_Normal.R
# https://github.com/alexanderrobitzsch/miceadds
# https://www.gerkovink.com/miceVignettes/Passive_Post_processing/Passive_imputation_post_processing.html


# Generate data -----------------------------------------------------------


baseline <- CauseSpecCovarMI::mds_shape_rates
shape_ev1 <- baseline[baseline$state == "REL", "shape"]
base_rate_ev1 <- baseline[baseline$state == "REL", "rate"]

ev1_pars <- list(
  "a1" = 1.5, 
  "h1_0" = 0.04,
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
  X_type = "continuous", 
  r = 0.5, 
  ev1_pars = ev1_pars,
  ev2_pars = ev2_pars, 
  rate_cens = baseline[baseline$state == "EFS", "rate"], 
  mech = "MAR", 
  p = 0.5,
  eta1 = -1
)

H1_nels <- dat$H1
H2_nels <- dat$H2
H1_true <- 0.04 * dat$t^(1.5)
H2_true <- baseline[baseline$state == "NRM", "rate"] * 
  dat$t^(baseline[baseline$state == "NRM", "shape"])

# Prepare imputations -----------------------------------------------------


# Set cumulative hazards and interactions as missing (they will be updated)
dat[, c("H1", "H2", "H1_Z", "H2_Z")] <- NA

# Make function to update basehaz
update_basehaz <- function(time, delta, x, z) {
  mod <- coxph(Surv(time, delta) ~ x + z, control = survival::coxph.control(timefix = FALSE))
  basehaz_df <- basehaz(mod, centered = FALSE)
  haz <- basehaz_df[match(time, basehaz_df[["time"]]), ][["hazard"]]
  return(haz)
}

# Prep matrices and methods
mat <- make.predictorMatrix(dat) 
mat[] <- 0L 
mat_ch1 <- mat_ch12 <- mat_ch12_int <- mat 
meth <- make.method(dat, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
meth_ch1 <- meth_ch12 <- meth_ch12_int <- meth 

# One cumhaz
mat_ch1["X", c("Z", "ev1", "H1")] <- 1
mat_ch1["H1", c("X", "Z", "t", "ev1", "ev2")] <- 1
meth_ch1["H1"] <- paste("~I(", expression(update_basehaz(t, ev1, X, Z)),")")
meth_ch1[c("H2", "H1_Z", "H2_Z")] <- ""

# Both cumhaz
mat_ch12["X", c("Z", "eps", "H1", "H2")] <- 1
mat_ch12["H1", c("X", "Z", "t", "ev1", "ev2")] <- 1
mat_ch12["H2", c("X", "Z", "t", "ev1", "ev2")] <- 1
meth_ch12["H1"] <- paste("~I(", expression(update_basehaz(t, ev1, X, Z)),")")
meth_ch12["H2"] <- paste("~I(", expression(update_basehaz(t, ev2, X, Z)),")")
meth_ch12[c("H1_Z", "H2_Z")] <- ""


# Imp settings
m <- c(15) # Number of imputations of interest
iters_MI <- 15

# Run imps
imp_ch1 <- mice::mice(
  dat, 
  m = m[length(m)],
  method = meth_ch1, 
  predictorMatrix = mat_ch1,
  maxit = iters_MI
)

imp_ch12 <- mice::mice(
  dat, 
  m = m[length(m)],
  method = meth_ch12, 
  predictorMatrix = mat_ch12,
  maxit = iters_MI
)

#make.method

imps_comp <- mice::complete(imp_ch12_int, action = "all")
cbind.data.frame(imps_comp$`1`$H1,
                 imps_comp$`2`$H1,
                 imps_comp$`3`$H1) |>  View()


# Make plots
plot(imp_ch12)
imps_comp <- mice::complete(imp_ch12, action = "all")


plot(
  x = H1_true, 
  y = H1_nels, 
  xlab = expression("True H"[10]*"(T)"),
  ylab = expression("Estimated H"[10]*"(T)"),
  xlim = c(0, 1.2), ylim = c(0, 1.2), type = "n"
)
abline(a = 0, b = 1, lwd = 2, col = "gray")
points(x = H1_true, y = H1_nels, col = "blue", pch = 16, cex = 0.5)
points(x = H1_true, y = imps_comp$`1`$H1, col = "black", pch = 16, cex = 0.5)
legend(
  "topleft", 
  legend = c("Iterative Breslow", "Nelson-Aalen"), 
  col = c("black", "blue"),
  pch = c(16, 16), 
  bty = "n"
)


plot(
  x = H2_true, 
  y = H2_nels, 
  xlab = expression("True H"[20]*"(T)"),
  ylab = expression("Estimated H"[20]*"(T)"),
  xlim = c(0, 1.2), ylim = c(0, 1.2), type = "n"
)
abline(a = 0, b = 1, lwd = 2, col = "gray")
points(x = H2_true, y = H2_nels, col = "blue", pch = 16, cex = 0.5)
points(x = H2_true, y = imps_comp$`1`$H2, col = "black", pch = 16, cex = 0.5)
legend(
  "topleft", 
  legend = c("Iterative Breslow", "Nelson-Aalen"), 
  col = c("black", "blue"),
  pch = c(16, 16), 
  bty = "n"
)

