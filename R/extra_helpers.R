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
  #' @importFrom ggplot2 ggplot aes geom_line geom_ribbon ylim theme
  #' scale_linetype_manual scale_color_manual ggtitle xlab ylab
  #' @importFrom tidyr replace_na
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
    mutate(times = .data$time)  %>% 
    left_join(CI_true, by = "times") %>% 
    
    # Correct for state at time zero
    replace_na(list(
      "true_pstate1" = 1, 
      "true_pstate2" = 0,
      "true_pstate3" = 0
    )) %>% 
    
    # Make into long format
    gather("state_est", "prob", .data$pstate1:.data$pstate3) %>% 
    gather("state_SE", "se", .data$se1:.data$se3) %>% 
    gather("state_true", "true", .data$true_pstate2:.data$true_pstate1) %>%
    unite(
      "state", 
      .data$state_est, 
      .data$state_SE, 
      .data$state_true
    ) %>% 
    mutate("state" = case_when(
      str_detect(.data$state, "pstate1_se1_true_pstate1") ~ "1",
      str_detect(.data$state, "pstate2_se2_true_pstate2") ~ "2",
      str_detect(.data$state, "pstate3_se3_true_pstate3") ~ "3"
    )) %>% 
    filter(!is.na(.data$state))%>%  
    
    # Add confidence intervals (these are plain)
    mutate(
      low = .data$prob - qnorm(0.975) * .data$se,
      low = ifelse(.data$low < 0, 0, .data$low),
      upp = .data$prob + qnorm(0.975) * .data$se,
      upp = ifelse(.data$upp > 1, 1, .data$upp)
    ) %>% 
    
    # Set up for labelling true and predicted
    gather("true_pred", "prob", .data$prob, .data$true) 
  
  
  # Plot begins here
  p1 <- ggplot(dat_plot,
               aes(.data$times, .data$prob, na.rm = T)) +
    geom_line(aes(col = .data$state, linetype = factor(.data$true_pred)), 
              size = 1.25) + ylim(c(0, 1)) +
    geom_ribbon(aes(x = .data$times, ymin = .data$low,
                    ymax = .data$upp, group = .data$state), 
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

