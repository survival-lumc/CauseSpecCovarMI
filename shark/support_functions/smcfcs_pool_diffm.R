########################################
## Pooling smcfcs() for a vector of m ##
########################################


# Function that:
# - Takes a smcfcs() object with m imputations eg. 100
# - Returns pooled estimates for that m, and other smaller m in the vector
# - Labels analysis type 
# - Depends on packages mitools, smcfcs
# - call source("quiet_printcat.R") if used on its own


smcfcs_pool_diffm <- function(imps, # output of smcfcs()
                              m, # vector of no. imputed datasets eg. c(1, 5, 10)
                              label,
                              true_betas) { # character vector labeling current analysis
  
  # Make list out of imputations
  impobj <- imputationList(imps$impDatasets)
  
  # Iterate procedure for all vals of m                     
  ests <- lapply(m, function(i) { # change indicator i later
    
    # Select first m imputed datasets
    if (i == 1) {
      impobj$imputations <- impobj$imputations[i]
      result <- with(impobj, coxph(Surv(t, eps == 1) ~ X1 + X2))[[1]]
      summ <- summary(result)$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")]
      
      # Format
      pooled <- cbind.data.frame(summ, confint(result)) %>%
        select(coef, se = "se(coef)", pval = "Pr(>|z|)", `2.5 %`, `97.5 %`)
      
    } else {
      impobj$imputations <- impobj$imputations[1:i]
      
      # Analyse and pool as usual 
      models <- with(impobj, coxph(Surv(t, eps == 1) ~ X1 + X2))
      summ <- quiet(summary(MIcombine(models)))
      
      # Summarise, add pva;
      pooled <-  summ %>% 
        mutate(z = results / se,
               pval = 2 * pnorm(abs(z), lower.tail = F)) %>%
        select(coef = results, se, pval, 
               `2.5 %` = `(lower`, `97.5 %` = `upper)`)
      
    }
    
    pooled <- pooled %>% 
      mutate(var = c("X1", "X2"), m = i, 
             analy = label, true = true_betas)
    
    return(pooled)
  })
  
  return(bind_rows(ests))
} 
