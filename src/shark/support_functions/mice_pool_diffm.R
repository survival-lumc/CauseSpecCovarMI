######################################
## Pooling mice() for a vector of m ##
######################################


# Function that:
# - Takes a mice object with m imputations eg. 100
# - Returns pooled estimates for that m, and other smaller m in the vector
# - Labels current replicate number and analysis type 
# - Depends on mice package


mice_pool_diffm <- function(imps, # output of mice()
                            m, # vector of no. imputed datasets eg. c(1, 5, 10)
                            j, # current replicate number
                            label,
                            true_betas) { # character vector labeling current analysis
  
  # Iterate procedure for all vals of m                     
  ests <- lapply(m, function(i) { # change indicator i later
    
    # Change the imp object so as to subsets first i imputations
    imps$m <- i
    
    # Have to use select() here, cannot use []
    imps$imp <- lapply(imps$imp, function(x) x %>% select(1:i))
    
    # Analyse and pool as usual 
    result <- with(imps, coxph(Surv(t, eps == 1) ~ X1 + X2))
    pooled <- suppressWarnings(summary(pool(result), conf.int = T)) 
    
    # If single imputation, selection is slightly different
    if (i == 1) {
      
      summ <- pooled$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")]
      pooled <- cbind.data.frame(summ, confint(result$analyses[[1]])) %>%
        select(coef, se = "se(coef)", pval = "Pr(>|z|)", `2.5 %`, `97.5 %`)
      
      # Single imputation so sd over imputations == 0
      #sd_imps <- rep(0, 2)
    } else {
      
      # Keep relevant estimates
      pooled <- as.data.frame(pooled) %>%
        select(coef = estimate, se = std.error, pval = p.value, `2.5 %`, `97.5 %`)
      
      # Compute sd of estimates across imputations
      #sd_imps <- apply(sapply(result$analyses, coefficients), 1, sd)
    }
    
    # Append all 
    pooled <- pooled %>%
      mutate(var = rownames(pooled), m = i, rep = j,
             analy = label, true = true_betas)
    
    return(pooled)
  })
  
  return(bind_rows(ests))
} 
