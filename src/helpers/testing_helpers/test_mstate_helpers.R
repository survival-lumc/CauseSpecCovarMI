##************************##
## Testing mstate helpers ##
##************************##

source("src/helpers/dat_generation_helpers.R")

# Work with complete data first, 

dat_MAR_bin <- generate_dat(n = 100,
                            X_type = "binary",
                            r = 0.5,
                            ev1_pars = list("a1" = 1, "h1_0" = 1, 
                                            "b1" = .5, "gamm1" = -.5),
                            ev2_pars = list("a2" = 1.7, "h2_0" = .5, 
                                            "b2" = -.5, "gamm2" =.5),
                            rate_cens = 1, 
                            mech = "MAR",
                            p = 0,
                            eta1 = 2)
                        

dat_MAR_contin <- generate_dat(n = 100,
                               X_type = "contin",
                               r = 0.5,
                               ev1_pars = list("a1" = 1, "h1_0" = 1, 
                                               "b1" = .5, "gamm1" = -.5),
                               ev2_pars = list("a2" = 1.7, "h2_0" = .5, 
                                               "b2" = -.5, "gamm2" =.5),
                               rate_cens = 1, 
                               mech = "MAR",
                               p = 0,
                               eta1 = 2)



library(mstate)

setup_mstate(dat_MAR_contin)

                    