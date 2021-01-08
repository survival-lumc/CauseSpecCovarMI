# Small sims to check whether interaction with of D with Z necessary

rep_num <- 1
base_seed <- 1

one_rep <- function(rep_num,
                    base_seed,
                    mpred) {
  
  set.seed(base_seed + rep_num)
  
  baseline <- readRDS(
    "./inst/testdata/MDS_shape_rates.rds"
  )
  
  shape_ev1 <- baseline[baseline$state == "REL", "shape"]
  base_rate_ev1 <- baseline[baseline$state == "REL", "rate"]
  
  
  # Parameter Weibull event 1
  ev1_pars <- list(
    "a1" = shape_ev1, 
    "h1_0" = base_rate_ev1,
    "b1" = c(0.75, 0.5), 
    "gamm1" = 1
  )
  
  # Parameters Weibull event 2
  ev2_pars <- list(
    "a2" = baseline[baseline$state == "NRM", "shape"], 
    "h2_0" = baseline[baseline$state == "NRM", "rate"], 
    "b2" = c(0.75, 1), 
    "gamm2" = .5
  )
  
  # Generate a dataset based on scenario
  dat <- SimsCauseSpecCovarMiss::generate_dat(
    n = 2000,
    X_type = "ordcat", 
    ev1_pars = ev1_pars,
    ev2_pars = ev2_pars, 
    rate_cens = baseline[baseline$state == "EFS", "rate"], 
    mech = "MAR", 
    p = 0.5,
    eta1 = -1
  )
  
  # Add interaction eps * Z
  inter_eps_Z <- as.data.frame(model.matrix(~ eps * Z, data = dat))
  dat$ev1_Z <- inter_eps_Z$`eps1:Z`
  dat$ev2_Z <- inter_eps_Z$`eps2:Z`
  
  # Make mice matrices 
  mpred_ch12 <- matrix(0, ncol(dat), ncol(dat),
                       dimnames = list(names(dat), 
                                       names(dat)))
  
  if (mpred == "inter") {
    mpred_ch12["X", c("Z", "eps", "H1", "H2", 
                      "ev1_Z", "ev2_Z")] <- 1
  } else {
    mpred_ch12["X", c("Z", "eps", "H1", "H2")] <- 1
  }
  
  
  # Get methods
  methods_mice <- SimsCauseSpecCovarMiss::set_mi_methods(
    dat = dat, 
    var_names_miss = naniar::miss_var_which(dat), 
    imp_type = "mice", 
  ) 
  
  
  imp_ch12 <- mice::mice(dat, m = 5,
                         method = methods_mice, 
                         predictorMatrix = mpred_ch12,
                         #maxit = 3, 
                         print = FALSE, 
                         threshold = 1) 
  
  res <- purrr::map(
    mice::complete(imp_ch12, action = "all"),
    ~ setup_mstate(.x)
  ) %$%
    summary(
      mice::pool(., dfcom = 999999), conf.int = T
    ) %>% 
    dplyr::select(-statistic, -df, -p.value)
    
  return(res)
}


library("furrr")
library(data.table)
library(magrittr)

plan(multisession, workers = 4)

test_nointer <- furrr::future_map_dfr(
  .x = 1:160,
  .f = ~ {
    devtools::load_all()
    one_rep(.x, base_seed = 209839, mpred = "nointer")
  },
  .id = "rep", 
  .progress = T
)

data.table(test_nointer) %>% 
  .[, .(est = mean(estimate)), by = term] %>% 
  print()



test_inter <- furrr::future_map_dfr(
  .x = 1:160,
  .f = ~ {
    devtools::load_all()
    one_rep(.x, base_seed = 10928, mpred = "inter")
  },
  .id = "rep", 
  .progress = T
)

data.table(test_inter) %>% 
  .[, .(est = mean(estimate)), by = term] %>% 
  print()