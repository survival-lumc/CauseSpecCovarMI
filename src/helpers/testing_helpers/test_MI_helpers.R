##**********************##
## Testing MICE helpers ##
##**********************##

source("src/helpers/dat_generation_helpers.R")
source("src/helpers/MI_helpers.R")


dat_MAR <- generate_dat(n = 500,
                        X_type = "contin",
                        r = 0.5,
                        ev1_pars = list("a1" = 1, "h1_0" = 1, 
                                        "b1" = .5, "gamm1" = -.5),
                        ev2_pars = list("a2" = 1.7, "h2_0" = .5, 
                                        "b2" = -.5, "gamm2" =.5),
                        rate_cens = 1, 
                        mech = "MCAR",
                        p = 0,
                        eta1 = 2)

dat_MAR$X_orig

mats <- get_predictor_mats(dat_MAR)
methods <- get_imp_models(dat_MAR)
m <- c(5, 10)

imp_ch1 <- mice(dat_MAR, m = m[length(m)],
                method = methods, 
                predictorMatrix = mats$CH12_int,
                maxit = 10)# 25
               #print = FALSE)

plot(imp_ch1)

complete(imp_ch1, 2)

lol <- complete(imp_ch1, action = "all")
hio <- lapply(lol, function(imp_dat) setup_mstate(imp_dat))

setup_mstate(dat_MAR %>% 
               mutate(X = X_orig))

# pipe with a dollar for lists!!
complete(imp_ch1, action = "all") %$%
  lapply(., function(imp_dat) setup_mstate(imp_dat)) %$%
  summary(pool(.))

summary(pool(hio))


