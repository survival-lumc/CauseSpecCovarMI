

generate_dat <- function(n,
                         X_type,
                         r) {
  
  #' @title Generates data.
  #' 
  #' @param n Sample size.
  #' @param X_type Either "binary", or "contin"
  #' @param r Desired correlation between X and Z.
  #' 
  #' @return Data-frame with missings induced
  
  # Computer true R needed if binary
  
  ifelse(X_type == "binary", 
         r <- pbiserial_to_pearson(p = 0.5, r_pb = r),
         r <- r)
  
  covmat <- matrix(c(1, r, 
                     r, 1), nrow = 2)
  
  dat_covars <- mvrnorm(n = n, mu = c(0, 0), Sigma = covmat) %>% 
    rename_all(c("X", "Z"))
  
  if (X_type == "binary") {
    
    dat_covars <- dat_covars %>% 
      mutate(X = ifelse(X <= x0, 0, 1))
    
  }
  
  return(dat_covars)
}



pbiserial_to_pearson <- function(p, r_pb) {
  
  #' @title Obtain necessary pearson to generate binary variable.
  #' 
  #' @param p Proportion of 1s in binary variable
  #' @param r_pb Desired point-biserial correlation
  #' 
  #' @return Pearson r to use when generating MVN data.
  
  x0 <- qnorm(p, lower.tail = F) # cutoff
  h <- dnorm(x0) # 'ordinate' = density/height of curve
  r <- r_pb * sqrt(p * (1 - p)) / h
  return(list("r" = r, "cutoff" = x0))
}
