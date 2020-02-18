##**********************************##
## Extra helpers (e.g. for visuals) ##
##**********************************##


cumincs_plot_truepred <- function(cox_long,
                                  combo,
                                  ev1_pars,
                                  ev2_pars,
                                  dat) {
  
  #' @title Visualise true/predicted probs with covariates
  #' 
  #' @inheritParams preds_mstate 
  #' @inheritParams get_true_cuminc
  #' @inheritParams get_predictor_mats
  #' 
  #' @return Ggplot by state and true/predicted
  #' 
  #' @export
  
  # Compute true
  CI_true <- get_true_cuminc(ev1_pars, ev2_pars, 
                             combo, times = dat$t)
  
  # Check vs probs trans
  new_dat <- data.frame(X.1 = c(combo$val_X, 0), 
                        X.2 = c(0, combo$val_X), 
                        Z.1 = c(combo$val_Z, 0), 
                        Z.2 = c(0, combo$val_Z), # THIS WAS SET TO VAL_X 
                        trans = c(1, 2), 
                        strata = c(1, 2))
  
  msfit_newdat <- msfit(cox_long, newdata = new_dat,
                        trans = trans.comprisk(2, c("Rel", "NRM")))
  
  # Make the plot
  dat_plot <- probtrans(msfit_newdat, predt = 0)[[1]] %>% 
    mutate(times = time)  %>% 
    left_join(CI_true, by = "times") %>% 
    
    # Correct for state at time zero
    replace_na(list(
      "true_pstate1" = 1, 
      "true_pstate2" = 0,
      "true_pstate3" = 0
    )) %>% 
    
    # Make into long format
    gather(state_est, prob, pstate1:pstate3) %>% 
    gather(state_SE, se, se1:se3) %>% 
    gather(state_true, true, true_pstate2:true_pstate1) %>%
    unite(state, state_est, state_SE, state_true) %>% 
    mutate(state = case_when(
      str_detect(state, "pstate1_se1_true_pstate1") ~ "1",
      str_detect(state, "pstate2_se2_true_pstate2") ~ "2",
      str_detect(state, "pstate3_se3_true_pstate3") ~ "3"
    )) %>% 
    filter(!is.na(state))%>%  
    
    # Add confidence intervals (these are plain)
    mutate(
      low = prob - qnorm(0.975) * se,
      low = ifelse(low < 0, 0, low),
      upp = prob + qnorm(0.975) * se,
      upp = ifelse(upp > 1, 1, upp)
    ) %>% 
    
    # Set up for labelling true and predicted
    gather(true_pred, prob, prob, true) 
  
  
  # Plot begins here
  p1 <- ggplot(dat_plot,
               aes(times, prob, na.rm = T)) +
    geom_line(aes(col = state, linetype = factor(true_pred)), 
              size = 1.25) + ylim(c(0, 1)) +
    geom_ribbon(aes(x = times, ymin = low, ymax = upp, group = state), 
                alpha = .25, col = NA) + 
    
    # Set up legends
    theme(legend.position = "bottom") +
    scale_linetype_manual(
      "Prob. type", 
      values = c("solid", "dotdash"),
      labels = c("Predicted", "True")
    ) +
    scale_color_manual(
      "State", 
      values = c(1, 2, 3),
      labels = c("EFS", "REL", "NRM")
    ) +
    
    # Title and labs
    ggtitle(paste0("Predicted and true probabilities for X = ",
                   combo[1],", Z = ", combo[2])) +
    xlab("Time") + 
    ylab("Probability")
  
  return(p1)
}

