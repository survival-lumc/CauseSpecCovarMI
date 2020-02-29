
dato <- mice::complete(imp_ch1, action = "all")[[1]]
mod_test <- setup_mstate(dato)


horiz <- c(1)

new_grid <- make_covar_grid(dat)[1, ]


# Would take 82 minutes per repetition with 13 covariate points
system.time(
hi <- preds_mstate(
  cox_long = mod_test,
  grid_obj = new_grid,
  times = horiz,
  ev1_pars = ev1_pars,
  ev2_pars = ev2_pars
)
)

tmat <- trans.comprisk(2, c("Rel", "NRM"))

new_dat <- data.frame(X.1 = c(new_grid$val_X, 0), 
                      X.2 = c(0, new_grid$val_X), 
                      Z.1 = c(new_grid$val_Z, 0), 
                      Z.2 = c(0, new_grid$val_Z), 
                      trans = c(1, 2), 
                      strata = c(1, 2))

msfit_newdat <- msfit(mod_test, newdata = new_dat,
                      trans = tmat)
View(msfit_newdat$Haz)


hio %>% 
  filter(time <= 1) %>% 
  slice(n()) %>% 
  select(time, pstate1, pstate2, pstate3) %>% 
  as.numeric()


lp_REL <- as.numeric(
  mod_test$coefficients[1:2] %*% c(new_grid$val_X, new_grid$val_Z)
)
lp_NRM <- as.numeric(
  mod_test$coefficients[3:4] %*% c(new_grid$val_X, new_grid$val_Z)
)


hio <- basehaz(mod_test, centered = F) %>% 
  pivot_wider(names_from = strata, 
              values_from = hazard) %>% 
  rename("H_REL" = `trans=1`,
         "H_NRM" = `trans=2`) %>% 
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

hio %>% 
  gather(state, prob, pstate1, pstate2, pstate3) %>% 
  ggplot(aes(time, prob, col = state)) +
  geom_line()



# try using survfit


estimates %>% 
  filter(m %in% c(0, 5)) %>% 
  ggplot(aes(analy, estimate)) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = true), linetype = "dashed") +
  geom_errorbar(aes(ymin = `2.5 %`,
                    ymax = `97.5 %`,
                    colour = analy), 
                width = 0.2,
                size = 1) +
  coord_flip() +
  facet_wrap(~ var) +
  theme(legend.position = "none")


