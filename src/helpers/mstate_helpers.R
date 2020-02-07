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


make_covar_grid <- function(.data) {
  
  sd_units <- c(0, 1, 2, -1, -2)
  z <- mean(.data$Z) + sd_units * sd(.data$Z)
  names(z) <- c("mean", "+1SD", "+2SD","-1SD","-2SD")
  
  if (is.factor(.data$X)) {
    
    x <- c(0, 1)
    names(x) <- c("0", "1")
    
    grid_obj <- expand.grid("X" = (names(x)),
                            "Z" = (names(z)), 
                            stringsAsFactors = F) %>% 
      as.data.frame() %>% 
      mutate(val_X = x[match(X, names(x))],
             val_Z = z[match(Z, names(z))])
                   
    # Maybe unite X and Z?
         
  } else {
    
    # Remove NA since there will be missings
    # check this?? Use X_orig instead?
    x <- mean(.data$X, na.rm = T) + sd_units * sd(.data$X, na.rm = T)
    names(x) <- names(z)
    
    grid_obj <- expand.grid("X" = names(x),
                            "Z" = names(z), 
                            stringsAsFactors = F) %>% 
      as.data.frame() %>% 
      unite("scen", c(X, Z), sep = "=") %>% 
      filter(!(str_detect(scen, "2SD") & 
                 !str_detect(scen, "mean"))) %>% 
      separate(scen, c("X", "Z"), sep = "=") %>%            
      mutate(val_X = x[match(X, names(x))],
             val_Z = z[match(Z, names(z))])
    
    # Maybe unite X and Z?
  }
  return(grid_obj)
}

make_covar_grid(dat_MAR) %>% 
  ggplot(aes(val_X, val_Z)) + geom_point()

library(mstate)
cox_long <- setup_mstate(dat_MAR)
grid_obj <- make_covar_grid(dat_MAR)
times <- c(0.2, 0.5)

# Do first on 1 grid obj
gridy <- grid_obj[1, ]


new_mds <- data.frame(X.1 = c(gridy$val_X, 0), # male
                      X.2 = c(0, gridy$val_X), # male
                      Z.1 = c(gridy$val_Z, 0), # 50 y.o
                      Z.2 = c(0, gridy$val_Z), # 50 y.o
                      trans = c(1, 2), 
                      strata = c(1, 2))


haz_weib <- Vectorize(function(alph, lam, t) {
  return(alph * lam * t^(alph - 1))
})

dens_weib <- Vectorize(function(alph, lam, t) {
  haz_weib(alph, lam, t) * exp(-lam * t^alph)
})

cumhaz_weib <- Vectorize(function(alph, lam, t) {
  lam * t^alph
})

# Take input vector of
gen_surv_weib <- Vectorize(function(cumhaz1, cumhaz2) {
  exp(-(cumhaz1 + cumhaz2))
})


# Make a general cumulative incidence function
cuminc_weib <- function(alph_ev, lam_ev, alph_comp, lam_comp, t) {
  
  prod <-  function(t) {
    haz_weib(alph_ev, lam_ev, t) * gen_surv_weib(cumhaz_weib(alph_ev, lam_ev, t),
                                                 cumhaz_weib(alph_comp, lam_comp, t))
  }
  
  ci_func <- Vectorize(function(upp) {
    integrate(prod, lower = 0, upper = upp)$value
  })
  
  return(ci_func(t))
}




preds_mstate <- function(cox_long,
                         grid_obj,
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
  
  # If not look at dat_gener_KMweib for true cuminc functions
  
  tmat <- trans.comprisk(2, c("Rel", "NRM"))
  
  predos <- lapply(1:nrow(grid_obj), function(row) {
    
    combo <- grid_obj[row, ]
    
    new_dat <- data.frame(X.1 = c(combo$val_X, 0), 
                          X.2 = c(0, combo$val_X), 
                          Z.1 = c(combo$val_Z, 0), 
                          Z.2 = c(0, combo$val_X), 
                          trans = c(1, 2), 
                          strata = c(1, 2))
    
    msfit_newdat <- msfit(cox_long, newdata = new_dat,
                          trans = tmat)
    
    preds <- probtrans(msfit_newdat, predt = 0)
    
    summ <- summary.probtrans(preds, times = times, conf.type = "log")[[1]]
    
    cbind.data.frame(summ, "X" = combo$X, "Z" = combo$Z)
  })
  
  return(bind_rows(predos))
}

preds_mstate(cox_long = setup_mstate(dat_MAR),
             grid_obj = make_covar_grid(dat_MAR), 
             times = c(0.2, 0.5))

cox_long <- setup_mstate(dat_MAR)
grid_obj <- 

# Work on expand grid thing here
#exp_grid_obj <- 



  
  

# Also need true cumulative incidences..