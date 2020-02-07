

missings <- miss_var_summary(dat_orig)

# Scenarios
beta <- c(0, 0.5, 1)
samp_size <- c(500, 2000)
miss <- c("low", "high")
X1_type <- c("contin", "binary")
miss_scen <- c("MCAR", "MAR_MRR", "MNAR",
               "MAR_gen")

scens <- expand.grid("beta" = beta, 
                     "miss" = miss, 
                     "X1_type" = X1_type, 
                     "scen" = miss_scen,
                     "n" = samp_size)

scens
