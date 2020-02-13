##************************##
## Testing mstate helpers ##
##************************##

devtools::load_all()

set.seed(123)

# Set weibull parameter globally as we will need them later
ev1_pars <- list("a1" = 2, "h1_0" = 1, 
                 "b1" = 1, "gamm1" = 1)

ev2_pars <- list("a2" = 2.5, "h2_0" = .5,
                 "b2" = .5, "gamm2" = .5)


# Lets start with a complete (large) dataset, with hardly any censoring
dat_full_contin <- generate_dat(n = 500,
                                X_type = "contin",
                                r = 0.1,
                                ev1_pars = ev1_pars,
                                ev2_pars = ev2_pars,
                                rate_cens = .01) 
                                #mech = "MNAR", 
                                #eta1 = 2, p = .5)

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
            
View(
  preds %>% 
    mutate_if(is.numeric, ~ round(., 3)) #%>% 
    #select(times, X, Z, 
     #      pstate1, pstate2, pstate3,
     #      true_pstate1, true_pstate2, true_pstate3)
)


# test correspondence!
combo <- 

cumincs_plot_truepred <- function(cox_long,
                                  combo,
                                  ev1_pars,
                                  ev2_pars,
                                  .data) {
  
  # Compute true
  CI_true <- get_true_cuminc(ev1_pars, ev2_pars, 
                             combo, times = .data$t)
  
  # Check vs probs trans
  new_dat <- data.frame(X.1 = c(combo$val_X, 0), 
                        X.2 = c(0, combo$val_X), 
                        Z.1 = c(combo$val_Z, 0), 
                        Z.2 = c(0, combo$val_Z), # THIS WAS SET TO VAL_X 
                        trans = c(1, 2), 
                        strata = c(1, 2))
  
  msfit_newdat <- msfit(cox_long, newdata = new_dat,
                        trans = trans.comprisk(2, c("Rel", "NRM")))
  
  prob_obj <- probtrans(msfit_newdat, predt = 0)[[1]] %>% 
    mutate(times = time) 
  
  probs <- prob_obj %>% 
    gather(state, prob, pstate1:pstate3) %>% 
    select(times, state, prob)
  
  preds <- prob_obj %>% 
    gather(state, se, se1:se3) %>% 
    mutate(state = case_when(
      str_detect(state, "1") ~ "pstate1",
      str_detect(state, "2") ~ "pstate2",
      str_detect(state, "3") ~ "pstate3"
    )) %>% 
    select(times, state, se) %>% 
    left_join(probs) %>% 
    mutate(
      low = prob - qnorm(0.975) * se,
      low = ifelse(low < 0, 0, low),
      upp = prob + qnorm(0.975) * se,
      upp = ifelse(upp > 1, 1, upp)
    )
  
  
  # Plot it
  p1 <- CI_true %>%
    gather(state, prob, true_pstate2:true_pstate1) %>% 
    bind_rows(preds) %>% 
    mutate(state_grp = case_when(
      str_detect(state, "1") ~ "1",
      str_detect(state, "2") ~ "2",
      str_detect(state, "3") ~ "3"
    )) %>% 
    mutate(true_pred = ifelse(str_detect(state, "true"), "true", "predicted")) %>% 
    group_by(state) %>% 
    fill(prob) %>% 
    ggplot(aes(times, prob)) +
    geom_line(aes(col = state_grp, linetype = factor(true_pred)), 
              size = 1.25) + ylim(c(0, 1)) +
    geom_ribbon(aes(x = times, ymin = low, ymax = upp, group = state), 
                alpha = .25, col = NA) + 
    theme_bw() +
    ggtitle(paste0("Predicted and true probabilities for X = ",
                   combo[1],", Z = ", combo[2])) +
    scale_linetype_manual("Prob. type", values = c("solid", "dotdash"),
                          labels = c("Predicted", "True")) +
    scale_color_manual("State", values = c(1, 2, 3),
                       labels = c("EFS", "REL", "NRM")) +
    theme(legend.position = "bottom", 
          plot.title = element_text(hjust = .5))
  
  return(p1)
}




cumincs_plot_truepred(cox_long = cox_long, 
                      combo = data.frame("val_X" = -2, "val_Z" = -2),
                      ev1_pars,
                      ev2_pars, 
                      dat_full_contin)















# Pipeline to apply preds_mstate to each imputed dataset ------------------






  
