##************************##
## Testing mstate helpers ##
##************************##

source("src/helpers/dat_generation_helpers.R")
source("src/helpers/mstate_helpers.R")
source("src/helpers/newsummaryprobtrans.r")

# Set weibull parameter globally as we will need them later
ev1_pars <- list("a1" = 2, "h1_0" = 1, 
                 #"b1" = .5, "gamm1" = 1)
                 "b1" = 0, "gamm1" = 0)

ev2_pars <- list("a2" = 2.5, "h2_0" = .5,
                 #"b2" = 1, "gamm2" = .5)
                 "b2" = 0, "gamm2" = 0)


# Lets start with a complete (large) dataset, with hardly any censoring
dat_full_contin <- generate_dat(n = 1000,
                                X_type = "contin",
                                r = 0,
                                ev1_pars = ev1_pars,
                                ev2_pars = ev2_pars,
                                rate_cens = .01, 
                                mech = "MAR",
                                p = 0,
                                eta1 = 2)

lapply(dat_full_contin, function(col) mean(is.na(col)))

table(dat_full_contin$eps)
hist(dat_full_contin$t)

# Run both cause-specific models using mstate
cox_long <- setup_mstate(dat_full_contin)
cox_long

# Make grid of covariate values at which to evaluate predictions
grid_obj <- make_covar_grid(dat_full_contin)
grid_obj

# Visualise it
grid_obj %>% 
  ggplot(aes(val_X, val_Z)) + geom_point()

# Time horizons
times <- c(0.5, 1)


preds <- preds_mstate(cox_long = cox_long,
                      grid_obj = grid_obj, 
                      times = times, 
                      ev1_pars = ev1_pars,
                      ev2_pars = ev2_pars)
            
view(
  preds %>% 
    mutate_if(is.numeric, ~ round(., 3)) %>% 
    select(times, X, Z, 
           pstate1, pstate2, pstate3,
           true_pstate1, true_pstate2, true_pstate3)
)



# Check empirically for marginal ones
combo <- data.frame("val_X" = 0, "val_Z" = 0)

CI_true <- get_true_cuminc(ev1_pars, ev2_pars, combo, times = dat_full_contin$t)

CI_emp <- survfit(Surv(t, eps) ~ 1, dat = dat_full_contin)

cbind.data.frame(CI_emp$time, CI_emp$pstate) %>% 
  rename_all(~ c("times", "pstate1", "pstate2", "pstate3")) %>% 
  left_join(CI_true) %>% 
  gather(state, prob, pstate1:true_pstate1) %>%  
  ggplot(aes(times, prob)) +
  geom_line(aes(col = state, linetype = state), size = 1.25) + ylim(c(0, 1)) +
  #geom_ribbon(aes(x = times, ymin = low, ymax = high, group = state), 
  #            alpha = .25, col = "grey") + 
  theme_bw()
  



# Check correspondence with covars and trans probs ------------------------




# Set weibull parameter globally as we will need them later
ev1_pars <- list("a1" = 2, "h1_0" = 1, 
                 "b1" = .5, "gamm1" = 1)
                 #"b1" = 0, "gamm1" = 0)

ev2_pars <- list("a2" = 2.5, "h2_0" = .5,
                 "b2" = 1, "gamm2" = .5)
                 #"b2" = 0, "gamm2" = 0)


# Lets start with a complete (large) dataset, with hardly any censoring
dat_full_contin <- generate_dat(n = 1000,
                                X_type = "contin",
                                r = 0,
                                ev1_pars = ev1_pars,
                                ev2_pars = ev2_pars,
                                rate_cens = .01, 
                                mech = "MAR",
                                p = 0,
                                eta1 = 2)

table(dat_full_contin$eps)
  
cox_long <- setup_mstate(dat_full_contin)

# Pick some covar combo
combo <- data.frame("val_X" = 2, "val_Z" = -2)


tmat <- trans.comprisk(2, c("Rel", "NRM"))

new_dat <- data.frame(X.1 = c(combo$val_X, 0), 
                      X.2 = c(0, combo$val_X), 
                      Z.1 = c(combo$val_Z, 0), 
                      Z.2 = c(0, combo$val_X), 
                      trans = c(1, 2), 
                      strata = c(1, 2))

msfit_newdat <- msfit(cox_long, newdata = new_dat,
                      trans = tmat)

preds <- probtrans(msfit_newdat, predt = 0)[[1]] %>% 
  mutate(times = time) %>% 
  select(times, pstate1, pstate2, pstate3)


# Compute true
CI_true <- get_true_cuminc(ev1_pars, ev2_pars, combo, times = dat_full_contin$t)

CI_true %>% 
  left_join(preds) %>% 
  fill(pstate1:pstate3) %>% 
  gather(state, prob, true_pstate2:pstate3) %>% 
  mutate(state_grp = case_when(
           str_detect(state, "1") ~ "1",
           str_detect(state, "2") ~ "2",
           str_detect(state, "3") ~ "3"
         )) %>% 
  ggplot(aes(times, prob)) +
  geom_line(aes(col = state_grp, linetype = state), 
            size = 1.25) + ylim(c(0, 1)) +
  #geom_ribbon(aes(x = times, ymin = low, ymax = high, group = state), 
  #            alpha = .25, col = "grey") + 
  theme_bw()
