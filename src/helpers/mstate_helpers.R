##************************##
## Testing mstate helpers ##
##************************##


setup_mstate <- function(.data) {
  
  #' @title Set-up mstate model (pre-predicting)
  #' 
  #' @description ...
  #' 
  #' @return Long cox model
  
  # Set up transition matrix 
  tmat <- trans.comprisk(2, c("Rel", "NRM"))
  covs <- c("X", "Z")
  
  # Long format
  dat_msprepped <- msprep(time = c(NA, "t", "t"),
                          status = c(NA, "ev1", "ev2"), 
                          data = .data,
                          trans = tmat,
                          keep = covs) 
  
  # Expand covariates
  dat_expanded <- expand.covs(dat_msprepped, covs,
                              append = TRUE, longnames = F)
  
  # Fit long cox model (both transitions)
  cox_long <- coxph(Surv(time, status) ~ 
                      X.1 + Z.1 + # Trans == 1
                      X.2 + Z.2 + # Trans == 2
                      strata(trans), # Separate baseline hazards
                    data = dat_expanded)
  
  # Prep objects for use with ms fit
  # Can you feed this into mice? Yes! summary(pool())
  
  return(cox_long)
}


preds_mstate <- function(cox_long,
                         exp_grid_obj,
                         times) {
  
  #' @title Predicting grids
  #' 
  #' @description ...
  #' 
  #' @param cox_long Cox model returned by setup_mstate
  #' @param times Prediction horizon(s) eg. c(2, 5, 10)
  #' 
  #' @return Long cox model
  
  #...
  
  
  
  new_mds <- data.frame(X.1 = c(0, 0), # male
                        X.2 = c(0, 0), # male
                        Z.1 = c(50, 0), # 50 y.o
                        Z.2 = c(0, 50), # 50 y.o
                        trans = c(1, 2), 
                        strata = c(1, 2))
  
  #data.frame(X.1 = c(1, 0), # male
  #           X.2 = c(0, 1), # male
  #           Z.1 = c(50, 0), # 50 y.o
  #          Z.2 = c(0, 50), # 50 y.o
  #           trans = c(1, 2), 
  #           strata = c(1, 2))

  
  mds_msfit <- msfit(cox_rel_long, newdata = new_mds,
                     trans = tmat)
  
  preds <- probtrans(mds_msfit, predt = 0)
  
  #summary.probtrans(preds, times = times, conf.type = "log")

  
  
}

# Work on expand grid thing here
#exp_grid_obj <- 


dat_MAR_contin

x <- mean(dat_MAR_contin$X) + c(0, 1, 2, -1, -2) * sd(dat_MAR_contin$X)
z <- mean(dat_MAR_contin$Z) + c(0, 1, 2, -1, -2) * sd(dat_MAR_contin$Z)

names(x) <- names(z) <- c("mean", "+1SD", "+2SD","-1SD","-2SD")

cums <- expand.grid("X" = names(x),
                    "Z" = names(z), 
                    stringsAsFactors = F)

grid_contin_names <- as.data.frame(cums) %>% 
  unite("scen", c(X, Z), sep = "=") %>% 
  filter(!(str_detect(scen, "2SD") & 
             !str_detect(scen, "mean"))) %>% 
  separate(scen, c("X", "Z"), sep = "=")

grid_contin_names %>% 
  mutate(X = x[match(X, names(x))],
         Z = z[match(Z, names(z))]) %>% 
  ggplot(aes(X, Z)) + geom_point()

grid_contin_names$x
match(grid_contin_names$x, names(x))

x[match(grid_contin_names$X, names(x))]

x <- c(0, 1)
z <- mean(dat_MAR_contin$Z) + c(0, 1, 2, -1, 2) * sd(dat_MAR_contin$Z)


# this is good for binary case!
names(x) <- x
names(z) <- c("mean", "+1SD", "+2SD","-1SD","-2SD")

cums <- expand.grid("x" = paste0("X_", names(x)),
                    "z" = paste0("Z_", names(z)), 
                    stringsAsFactors = F)

cums

paste0("X_", names(x))

# Also need true cumulative incidences..