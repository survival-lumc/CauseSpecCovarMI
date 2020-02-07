##**************************##
## Testing dat_generation.R ##
##**************************##


# Visualise with continous X1, change below to "binary"
# for equivalent visualisations
# Note if shape of one of two events is too small then
# marginal cum base haz will have some NAs at the beginning
# due to time points not being 'distinct' enough (eg 10e-9 and 10e-8)

source("src/helpers/dat_generation_helpers.R")


# MCAR --------------------------------------------------------------------


dat_MCAR <- generate_dat(n = 1000, # 1000
                         X_type = "contin",
                         r = 0.5,
                         ev1_pars = list("a1" = 1, "h1_0" = 1, 
                                         "b1" = .5, "gamm1" = -.5),
                         ev2_pars = list("a2" = 1.7, "h2_0" = .5, 
                                         "b2" = -.5, "gamm2" =.5),
                         rate_cens = 1, 
                         mech = "MCAR",
                         p = 0.5)

# Check correlation and missingness proportion
mean(is.na(dat_MCAR$X))
cor(dat_MCAR$X_orig, dat_MCAR$Z)
    
# Check cause-specific hazards models
# library(survival)
coxph(Surv(t, eps == 1) ~ X_orig + Z, data = dat_MCAR)
coxph(Surv(t, eps == 2) ~ X_orig + Z, data = dat_MCAR)

# Visualise missingness
ggplot(dat_MCAR, aes(X_orig, Z, col = factor(miss_ind))) +
  geom_point() +
  scale_colour_manual("Missing indicator",
                      labels = c("Observed", "Missing"),
                      values = c(9, 6))


# MAR (depending on Z) ----------------------------------------------------


dat_MAR <- generate_dat(n = 1000,
                        X_type = "contin",
                        r = .7,
                        ev1_pars = list("a1" = 1, "h1_0" = 1, 
                                        "b1" = .5, "gamm1" = -.5),
                        ev2_pars = list("a2" = 1.7, "h2_0" = .5, 
                                        "b2" = -.5, "gamm2" =.5),
                        rate_cens = 1, 
                        mech = "MAR",
                        p = 0.5,
                        eta1 = -.5)

ggplot(dat_MAR, aes(X_orig, Z, col = factor(miss_ind))) +
  geom_point() +
  scale_colour_manual("Missing indicator",
                      labels = c("Observed", "Missing"),
                      values = c(9, 6))
                    


# MNAR (depending on X itself) --------------------------------------------


dat_MNAR <- generate_dat(n = 1000,
                         X_type = "contin",
                         r = 0.5,
                         ev1_pars = list("a1" = 1, "h1_0" = 1, 
                                         "b1" = .5, "gamm1" = -.5),
                         ev2_pars = list("a2" = 1.7, "h2_0" = .5, 
                                         "b2" = -.5, "gamm2" =.5),
                         rate_cens = 1, 
                         mech = "MNAR",
                         p = 0.5,
                         eta1 = 1)
                        
ggplot(dat_MNAR, aes(X_orig, Z, col = factor(miss_ind))) +
  geom_point() +
  scale_colour_manual("Missing indicator",
                      labels = c("Observed", "Missing"),
                      values = c(9, 6))


# MAR_GEN (depending on t, long follow has less missing) ------------------



dat_MAR_GEN <- generate_dat(n = 1000,
                            X_type = "contin",
                            r = 0.5,
                            ev1_pars = list("a1" = 1, "h1_0" = 1, 
                                            "b1" = .5, "gamm1" = 1),
                            ev2_pars = list("a2" = 1.7, "h2_0" = .5, 
                                            "b2" = -.5, "gamm2" = 1),
                            rate_cens = 1, 
                            mech = "MAR_GEN",
                            p = 0.5,
                            eta1 = -2)
                         

ggplot(dat_MAR_GEN, aes(X_orig, Z, col = factor(miss_ind))) +
  geom_point() +
  #ylab("log(t)") +
  scale_colour_manual("Missing indicator",
                      labels = c("Observed", "Missing"),
                      values = c(9, 6))

library(naniar)

# So earlier t more likely to be missing
ggplot(dat_MAR_GEN, aes(t, X)) +
  geom_miss_point()




# Binary, point biserial and logreg ---------------------------------------



dat_bin <- generate_dat(n = 1000, # 1000
                        X_type = "binary",
                        r = 0.75,
                        ev1_pars = list("a1" = 1, "h1_0" = 1, 
                                        "b1" = .5, "gamm1" = -.5),
                        ev2_pars = list("a2" = 1.7, "h2_0" = .5, 
                                        "b2" = -.5, "gamm2" =.5),
                        rate_cens = 1, 
                        mech = "MCAR",
                        p = 0)
                         

# Check correlation and missingness proportion
mean(is.na(dat_bin$X))
cor(as.numeric(dat_bin$X_orig), dat_bin$Z)

linmod <- glm(as.numeric(X_orig) ~ Z, data = dat_bin)

logreg <- glm(X_orig ~ Z, data = dat_bin, family = binomial)

plot(dat_bin$Z, dat_bin$X_orig)
lines(dat_bin$Z, predict(linmod))

plot(dat_bin$Z, plogis(predict(logreg)), ylim = c(-1, 2))
points(dat_bin$Z, as.numeric(dat_bin$X_orig) - 1)

