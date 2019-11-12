

format_MI_repet <- function(repetitions) {
  
  # Extract coef and se of first repetition
  se_rep1 <- bind_rows(repetitions) %>% 
    filter(rep == 1) %>% 
    select(coef, se, var, analy, m) %>% 
    rename(coef_i1 = coef)
  
  # Bring the rest together
  agg <- bind_rows(repetitions) %>%
    
    # Power and coverage
    mutate(pow = pval < 0.05,
           cover = `2.5 %` < true & true < `97.5 %`,
           bias = coef - true) %>%
    group_by(var, analy, m) %>%
    mutate(sd_reps = sd(coef)) %>%
    summarise_all(~ round(mean(.), 3)) %>%
    ungroup() %>%
    select(-rep, -se) %>%
    
    # Join se
    right_join(se_rep1, by = c("var", "analy", "m")) %>%
    as.data.frame()
  
  return(agg)
}


format_CCA_ref <- function(mod, true_betas, label) {
  
  summ_mod <- summary(mod)$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")]
  results <- cbind.data.frame(summ_mod, confint(mod)) %>%
    rownames_to_column("var") %>% 
    mutate(var = ifelse(var == "X1_orig", "X1", var),
           analy = "label",
           true = true_betas,
           m = 0, pval = `Pr(>|z|)`,
           sd_reps = 0,
           pow = pval < 0.05,
           cover = `2.5 %` < true & true < `97.5 %`,
           bias = coef - true,
           coef_i1 = coef) %>% 
    rename(se = `se(coef)`)

  return(results)
}
