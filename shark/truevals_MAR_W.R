###########################################
## Get reference values for MAR with W   ##
## by simulating large data (n = 100000) ##
###########################################


# For reproducibility
set.seed(1984)

# For beta1 = 0
dat_b0 <- dat_gener_W(N = 100000, 
                      X_type = "contin",
                      mus = c(0, 0), 
                      covmat = matrix(c(1, 0.25, 
                                        0.25, 1), nrow = 2), # make into correlation mat
                      mech = "MCAR", 
                      pars_MAR = c(1, 1),
                      p = 0.1, 
                      cause2 = "weib", 
                      vals_t1 = c(0.3, 1, c(0, -0.5)), 
                      vals_t2 = c(1.7, 0.5, -0.5, 0.5)) 


cox_b0 <- coxph(Surv(t, eps == 1) ~ X1_orig + X2, data = dat_b0)
cox_b0$coefficients


# For beta1 = 0.5
dat_b0.5 <- dat_gener_W(N = 100000, 
                        X_type = "contin",
                        mus = c(0, 0), 
                        covmat = matrix(c(1, 0.25, 
                                          0.25, 1), nrow = 2), # make into correlation mat
                        mech = "MCAR", 
                        pars_MAR = c(1, 1),
                        p = 0.1, 
                        cause2 = "weib", 
                        vals_t1 = c(0.3, 1, c(0.5, -0.5)), 
                        vals_t2 = c(1.7, 0.5, -0.5, 0.5)) 


cox_b0.5 <- coxph(Surv(t, eps == 1) ~ X1_orig + X2, data = dat_b0.5)
cox_b0.5$coefficients


# For beta1 = 1
dat_b1 <- dat_gener_W(N = 100000, 
                      X_type = "contin",
                      mus = c(0, 0), 
                      covmat = matrix(c(1, 0.25, 
                                        0.25, 1), nrow = 2), # make into correlation mat
                      mech = "MCAR", 
                      pars_MAR = c(1, 1),
                      p = 0.1, 
                      cause2 = "weib", 
                      vals_t1 = c(0.3, 1, c(1, -0.5)), 
                      vals_t2 = c(1.7, 0.5, -0.5, 0.5)) 


cox_b1 <- coxph(Surv(t, eps == 1) ~ X1_orig + X2, data = dat_b1)
cox_b1$coefficients
