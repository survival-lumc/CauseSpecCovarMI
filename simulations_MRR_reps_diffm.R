# Attempt with different m and data sets 

source("dat_generation.R")

set.seed(1984)

N <- 5 # number of datasets to create
true_betas <- c(0.5, -0.5) #beta and gamma
mech <- "MAR"
m <- c(1, 5, 10)
method <- "norm"
prob <- 0.5


# Generate independent datasets & assign index
sims <- lapply(1:N, function(i) {
  
  dat <- dat_gener(N = 200, 
                   X_type = "contin",
                   mus = c(0, 0), 
                   covmat = matrix(c(1, 0.25, 
                                     0.25, 1), nrow = 2), 
                   mech = mech, 
                   p = prob, 
                   cause2 = "weib", 
                   vals_t1 = c(0.3, 1, true_betas), 
                   vals_t2 = c(1.7, 0.5, -0.5, 0.5)) 
  
  return(list("dat" = dat, "rep" = i))
})



R <- 10 # umber p

# 50 secs per dataset, so for 1000 thats 14h (so m = c(1, 10, 50) and R = 10)
system.time(
final_test <- lapply(sims, function(obj) {
  
  # Extract data set
  dat <- obj$dat
  
  # Set seed
  set.seed(obj$rep)
  
  testips <- lapply(1:R, function(r) {
    
    set.seed(r)
    
    ## Multiple imputation section:
    
    # Make predictor matrices for both ch1 and ch12
    mpred_ch1 <- mpred_ch12 <- mpred_ch12_int <-  matrix(0, ncol(dat), ncol(dat),
                                                         dimnames = list(names(dat), names(dat)))
    
    
    ## Ch1 model: MI with Z, eps (as factor) and H1(t) as covars
    mpred_ch1["X", c("Z", "ev1", "H1")] <- 1
    
    
    # Run m = 10 imputations, run cause-specific cox, and pool - 
    # tidyversify this using https://stefvanbuuren.name/fimd/workflow.html
    imp_ch1 <- mice(dat, m = m[length(m)],
                    method = method, 
                    predictorMatrix = mpred_ch1,
                    print = FALSE)
    
    results_ch1 <- bind_rows(single_rep_m(imps = imp_ch1,
                                          m = m, r = r, label = "ch1"))
    
    
    ## Ch12 model: MI with Z, eps, H1(t) & H2(t) as covars
    mpred_ch12["X", c("Z", "eps", "H1", "H2")] <- 1
    
    # Run m = 10 imputations, run cause-specific cox, and pool
    imp_ch12 <- mice(dat, m = m[length(m)],
                     method = method, 
                     predictorMatrix = mpred_ch12,
                     print = FALSE)
    
    results_ch12 <- bind_rows(single_rep_m(imps = imp_ch12,
                                           m = m, r = r, label = "ch12"))
    
    
    # Ch12 model with addition of interactions H x Z
    mpred_ch12_int["X", c("Z", "eps", "H1", "H2", "H1_Z", "H2_Z")] <- 1
    
    imp_ch12_int <- mice(dat, m = m[length(m)],
                         method = method,
                         predictorMatrix = mpred_ch12_int,
                         print = FALSE)
    
    results_ch12_int <- bind_rows(single_rep_m(imps = imp_ch12_int, 
                                               m = m, r = r, label = "ch12_int"))
    
    results_MI <- rbind.data.frame(results_ch1, results_ch12, results_ch12_int)
    return(results_MI) 
  })
  
  
  # tapply with multiple indices? Instead of group_by()
  # Also how to do coverage/power? 
  agg <- bind_rows(testips) %>%
    group_by(var, analy, m) %>%
    mutate(sd_reps = sd(coef)) %>%
    summarise_all(~ round(mean(.), 3)) %>%
    ungroup() %>%
    select(-rep) %>%
    as.data.frame()
  
  
  ## Ref model: 
  mod_ref <- coxph(Surv(t, eps == 1) ~ X_orig + Z_orig, data = dat)
  
  # Read in summary, compute CI & label analysis type
  summ_ref <- summary(mod_ref)$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")]
  results_ref <- cbind.data.frame(summ_ref, confint(mod_ref)) %>%
    mutate(analy = "1ref")
  
  
  ## CCA:
  mod_CCA <- coxph(Surv(t, eps == 1) ~ X + Z, data = dat)
  
  # Read in summary
  summ_CCA <- summary(mod_CCA)$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")]
  
  # Combine / format results together with reference
  results_CCA <- cbind.data.frame(summ_CCA, confint(mod_CCA)) %>%
    mutate(analy = "CCA") %>% 
    bind_rows(results_ref) %>%
    mutate(var = rep(c("X", "Z"), 2),
           true = rep(true_betas, 2),
           m = 0, pval = `Pr(>|z|)`,
           sd_imps = 0,
           sd_reps = 0) %>%
    select(var, analy, m, coef, se = `se(coef)`, 
           pval, `2.5 %`, `97.5 %`, sd_imps, true, sd_reps) 
  
  return(rbind.data.frame(results_CCA, agg))
})
)

m <- c(1, 5, 10)
N <- 1

bind_rows(final_test) %>%
  group_by(var, analy, m) %>%
  mutate(se_emp = sd(coef)) %>%
  summarise_all(~ round(mean(.), 3)) %>%
  select(var, analy, m, coef, se, se_emp) %>%
  as.data.frame()




single_rep_m <- function(imps, # output of mice()
                         m, # vector of no. imputed datasets eg. c(1, 5, 10)
                         r, # current replicate number
                         label) { # character vector labeling current analysis
  
  # Iterate procedure for all vals of m                     
  ests <- lapply(m, function(i) { # change indicator i later
    
    # Change the imp object so as to subsets first i imputations
    imps$m <- i
    
    # Have to use select() here, cannot use []
    imps$imp <- lapply(imps$imp, function(x) x %>% select(1:i))
    
    # Analyse and pool as usual 
    result <- with(imps, coxph(Surv(t, eps == 1) ~ X + Z))
    pooled <- suppressWarnings(summary(pool(result), conf.int = T)) 
    
    # If single imputation, selection is slightly different
    if (i == 1) {
      
      summ <- pooled$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")]
      pooled <- cbind.data.frame(summ, confint(result$analyses[[1]])) %>%
        select(coef, se = "se(coef)", pval = "Pr(>|z|)", `2.5 %`, `97.5 %`)
      
      # Single imputation so sd over imputations == 0
      sd_imps <- rep(0, 2)
    } else {
      
      # Keep relevant estimates
      pooled <- as.data.frame(pooled) %>%
        select(coef = estimate, se = std.error, pval = p.value, `2.5 %`, `97.5 %`)
      
      # Compute sd of estimates across imputations
      sd_imps <- apply(sapply(result$analyses, coefficients), 1, sd)
    }
    
    # Append all 
    pooled <- pooled %>%
      mutate(var = c("X", "Z"), m = i, sd_imps, rep = r,
             analy = label, true = true_betas)
  
    return(pooled)
  })
} 
