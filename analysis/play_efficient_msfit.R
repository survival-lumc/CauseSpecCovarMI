################################################
## Efficiently set up predicted probabilities ##
################################################


# We only set up fo
preds_mstate <- function(cox_long,
                         grid_obj,
                         times,
                         ev1_pars,
                         ev2_pars) {
  
  #' @title Predicting grids
  #' 
  #' @description ...
  #' 
  #' @param cox_long Cox model returned by setup_mstate
  #' @param times Prediction horizon(s) eg. c(2, 5, 10)
  #' @inheritParams get_true_cuminc
  #' @param grid_obj Grid object
  #' 
  #' @importFrom mstate msfit probtrans
  #' @importFrom purrr map modify_depth map_dfr imap_dfr
  #' 
  #' @return Long cox model
  #' 
  #' @export
  
  #...
  
  # If not look at dat_gener_KMweib for true cuminc functions
  
  tmat <- trans.comprisk(2, c("Rel", "NRM"))
  
  purrr::map_dfr(1:nrow(grid_obj), function(row) {
    
    combo <- grid_obj[row, ]
    
    new_dat <- data.frame(X.1 = c(combo$val_X, 0), 
                          X.2 = c(0, combo$val_X), 
                          Z.1 = c(combo$val_Z, 0), 
                          Z.2 = c(0, combo$val_Z), 
                          trans = c(1, 2), 
                          strata = c(1, 2))
    
    msfit_newdat <- msfit(cox_long, newdata = new_dat,
                          trans = tmat)
    
    preds <- record_warning(probtrans(msfit_newdat, predt = 0))
    
    
    # Add log later after pooling
    summ <- summary.probtrans(preds$value, 
                              times = times, conf.type = "none")[[1]]
    
    # Add the true ones
    true_CI <- get_true_cuminc(ev1_pars = ev1_pars, 
                               ev2_pars = ev2_pars,
                               combo = combo, 
                               times = times)
    
    cbind.data.frame(summ, 
                     "X" = combo$X, 
                     "Z" = combo$Z) %>% 
      left_join(true_CI, by = "times") %>% 
      
      # Record 1 if probtrans error was recorded
      mutate(
        warning = ifelse(preds$warning != "0", 1, 0)
      ) %>% 
      dplyr::mutate_if(is.factor, as.character) # avoid bind_rows warnings
  })
}






dato <- mice::complete(imp_ch12_int, action = "all")[[2]]
mod_test <- setup_mstate(dato)

new_dat <- data.frame(X.1 = c(0, 0), 
                      X.2 = c(0, 0), 
                      Z.1 = c(0, 0), 
                      Z.2 = c(0, 0), 
                      trans = c(1, 2), 
                      strata = c(1, 2))

# Compute baseline
tmat <- trans.comprisk(2, c("Rel", "NRM"))
msfit_newdat <- msfit(mod_test, 
                      newdata = new_dat,
                      trans = tmat)

# This is covariate combo to evaluation
new_grid <- make_covar_grid(dat)[12, ]

new_grid

# Compute linear predictors
lp_REL <- as.numeric(
  mod_test$coefficients[1:2] %*% c(new_grid$val_X, new_grid$val_Z)
)

lp_NRM <- as.numeric(
  mod_test$coefficients[3:4] %*% c(new_grid$val_X, new_grid$val_Z)
)


# Use 
preds_test <- msfit_newdat$Haz %>% 
  pivot_wider(names_from = trans, 
              values_from = Haz) %>% 
  rename("H_REL" = `1`,
         "H_NRM" = `2`) %>% 
  mutate(H_REL = H_REL * exp(lp_REL),
         H_NRM = H_NRM * exp(lp_NRM)) %>% 
  mutate(pstate1 = exp(-(H_REL + H_NRM))) %>% 
  add_row(time = 0, H_REL = 0, H_NRM = 0, pstate1 = 1) %>% 
  arrange(time) %>% 
  mutate(haz_REL = c(0, diff(H_REL)),
         haz_NRM = c(0, diff(H_NRM)),
         EFS_min1 = c(0, pstate1[-length(pstate1)]),
         pstate2 = cumsum(haz_REL * EFS_min1),  
         pstate3 = cumsum(haz_NRM * EFS_min1),
         check = pstate1 + pstate2 + pstate3) 

horiz <- c(1, 5, 10)

map_dfr(horiz, 
        ~ preds_test %>% 
          filter(time <= .x) %>% 
          slice(n()) %>% 
          select(time, pstate1, pstate2, pstate3) %>% 
          as.data.frame())


# Compare results with probtrans ------------------------------------------


new_dat_probtrans <- data.frame(X.1 = c(new_grid$val_X, 0),
                                X.2 = c(0, new_grid$val_X), 
                                Z.1 = c(new_grid$val_Z, 0), 
                                Z.2 = c(0, new_grid$val_Z), 
                                trans = c(1, 2), 
                                strata = c(1, 2))
                      


msfit_newdat_probtrans <- msfit(mod_test, 
                                newdata = new_dat_probtrans,
                                trans = tmat)

probtrans_test <- probtrans(msfit_newdat_probtrans, predt = 0)

summary.probtrans(probtrans_test, times = 1, conf.type = "none")[[1]] %>% 
  as.data.frame() %>% 
  select(times, pstate1, pstate2, pstate3) %>% 
  round(., 4)

preds_test %>% 
  filter(time <= 1) %>% 
  slice(n()) %>% 
  select(time, pstate1, pstate2, pstate3) %>% 
  as.data.frame() %>% 
  round(., 4)


preds_mstate2 <- function(cox_long,
                          grid_obj,
                          times,
                          ev1_pars,
                          ev2_pars) {
  
  # If not look at dat_gener_KMweib for true cuminc functions
  
  tmat <- trans.comprisk(2, c("Rel", "NRM"))
  
  new_dat_probtrans <- data.frame(
    X.1 = c(0, 0),
    X.2 = c(0, 0), 
    Z.1 = c(0, 0), 
    Z.2 = c(0, 0), 
    trans = c(1, 2), 
    strata = c(1, 2)
  )
  
  msfit_newdat_probtrans <- msfit(
    mod_test, 
    newdata = new_dat_probtrans,
    trans = tmat
  )
  
  # Do it for all grid points
  purrr::map_dfr(1:nrow(grid_obj), function(row) {
    
    combo <- grid_obj[row, ]
    
    # Compute linear predictors
    lp_REL <- as.numeric(
      cox_long$coefficients[1:2] %*% c(combo$val_X, combo$val_Z)
    )
    
    lp_NRM <- as.numeric(
      cox_long$coefficients[3:4] %*% c(combo$val_X, combo$val_Z)
    )
    
    # This is where the magic happens
    preds_test <- msfit_newdat$Haz %>% 
      pivot_wider(names_from = trans, 
                  values_from = Haz) %>% 
      rename("H_REL" = `1`,
             "H_NRM" = `2`) %>% 
      mutate(H_REL = H_REL * exp(lp_REL),
             H_NRM = H_NRM * exp(lp_NRM)) %>% 
      mutate(pstate1 = exp(-(H_REL + H_NRM))) %>% 
      add_row(time = 0, H_REL = 0, H_NRM = 0, pstate1 = 1) %>% 
      arrange(time) %>% 
      mutate(haz_REL = c(0, diff(H_REL)),
             haz_NRM = c(0, diff(H_NRM)),
             EFS_min1 = c(0, pstate1[-length(pstate1)]),
             pstate2 = cumsum(haz_REL * EFS_min1),  
             pstate3 = cumsum(haz_NRM * EFS_min1),
             check = pstate1 + pstate2 + pstate3) 
    
    
    # Find predicions at time horizons
    summ <- map_dfr(times, 
                    ~ preds_test %>% 
                      filter(time <= .x) %>% 
                      slice(n()) %>% 
                      select(times = time, 
                             pstate1, pstate2, pstate3) %>% 
                      mutate(times = ceiling(times)) %>% 
                      as.data.frame())
    
    # Add the true ones
    true_CI <- get_true_cuminc(ev1_pars = ev1_pars, 
                               ev2_pars = ev2_pars,
                               combo = combo, 
                               times = times)
    
    cbind.data.frame(summ, 
                     "X" = combo$X, 
                     "Z" = combo$Z) %>% 
      left_join(true_CI, by = "times") %>% 
      dplyr::mutate_if(is.factor, as.character) # avoid bind_rows warnings
  })
}


preds_mstate2(mod_test, 
              make_covar_grid(dato), 
              times = c(1, 5, 10), 
              ev1_pars = ev1_pars, 
              ev2_pars = ev2_pars) %>% 
  select(times, 
         pstate1, true_pstate1, 
         pstate2, true_pstate2, 
         pstate3, true_pstate3)



dat %>% 
  gather(ev, cumhaz, H1:H2) %>% 
  ggplot(aes(t, cumhaz, col = ev)) +
  geom_line(size = 1.5)

dat %>% 
  ggplot(aes(H1, H2)) + geom_point()

