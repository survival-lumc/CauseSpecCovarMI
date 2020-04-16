##************************************##
##   Importing simulation results     ##
##************************************##


# Load compendium 
devtools::load_all()


# Import estimates --------------------------------------------------------


# * Process pilots --------------------------------------------------------


# Make pattern corresponding to 14 pilot scenarios
# This is necessary since they were not in same format
pilot_nums <- paste(
  sapply(as.character(1:14), function(num) paste0("scen", num, "_")),
  collapse = "|"
)

# Keep labels of scenario for separating "scen_summary" colums
labs_scens <- c("n", "prop_miss", "beta1", "miss_mech", "X_level",
                "rho", "eta1", "haz_shape")

# Import and format
pilot_estims <- list.files(
  path = "./analysis/simulations/sim-reps_summarised/estimates/",
  pattern = paste0("(", pilot_nums, ")"), 
  full.names = T
) %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist() %>% 
  .[, ':=' (
    rmse = NULL,
    mcarlo_se_cover = sqrt((cover * (1 - cover)) / n),
    scen_summary = paste0(scen_summary, "-haz_shape=similar")
  )]


# * Process and bind rest -------------------------------------------------


# Make list of all files 
all_files <- list.files(
  path = "./analysis/simulations/sim-reps_summarised/estimates/",
  full.names = T
)

# Exclude pilots using grep
all_estims <- grep(
  x = all_files,
  pattern = paste0("(", pilot_nums, ")"), 
  invert = T, 
  value = T
)  %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist() %>% 
  
  # Bind pilots
  rbind(., pilot_estims) %>% 
  
  # Adjust for eta minus label
  .[, scen_summary := gsub(
    "eta1=-", 
    "eta1=min", 
    scen_summary
  )] %>%

  # Separate the scenario columns
  .[, (labs_scens) := data.table::tstrsplit(scen_summary, split = "-")] %>% 

  # Relevel methods
  .[, analy := factor(
    analy, levels = c("ref", "CCA", "ch1", "ch12", "ch12_int", "smcfcs")
  )] %>% 
  
  # Relevel other factors and properly names
  .[, ':=' (
    prop_miss = factor(
      prop_miss, 
      levels = c("prop_miss=0.1","prop_miss=0.5"),
      labels = c("10%", "50%")
    ),
    haz_shape = factor(
      haz_shape,
      levels = c("haz_shape=similar", "haz_shape=different"),
      labels = c("similar", "different")
    ),
    beta1 = factor(
      beta1,
      levels = c("beta1=0", "beta1=0.5", "beta1=1"),
      labels = c("0", "0.5", "1")
    ),
    eta1 = factor(
      eta1,
      levels = c("eta1=NA", "eta1=min1", "eta1=min2"),
      labels = c("None", "Weak", "Strong")
    ),
    miss_mech = factor(
      miss_mech,
      levels = c("miss_mech=MCAR", "miss_mech=MAR", 
                 "miss_mech=MAR_GEN", "miss_mech=MNAR"),
      labels = c("MCAR", "MAR", "MAR_GEN", "MNAR")
    ),
    X_level = factor(
      X_level,
      levels = c("X_level=continous", "X_level=binary"),
      labels = c("continuous", "binary")
    ),
    rho = factor(
      rho,
      levels = "rho=0.5",
      labels = "0.5"
    ),
    n = factor(
      n,
      levels = "n=2000",
      labels = "2000"
    )
  )]

rm(pilot_estims, all_files)


# Import predictions ------------------------------------------------------


# * Process pilots --------------------------------------------------------


# Import and format
pilot_preds <- list.files(
  path = "./analysis/simulations/sim-reps_summarised/predictions/",
  pattern = paste0("(", pilot_nums, ")"), 
  full.names = T
) %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist() %>% 
  .[, scen_summary := paste0(scen_summary, "-haz_shape=similar")]


# * Process and bind rest -------------------------------------------------


# Make list of all files 
all_files <- list.files(
  path = "./analysis/simulations/sim-reps_summarised/predictions/",
  full.names = T
)

# Exclude pilots using grep
all_preds <- grep(
  x = all_files,
  pattern = paste0("(", pilot_nums, ")"), 
  invert = T, 
  value = T
)  %>% 
  purrr::map(readRDS) %>% 
  data.table::rbindlist() %>% 

  # Exclude rmse_se, and mcse_rmse (for now)
  .[, ':=' (
    rmse_se = NULL,
    mcarlo_se_rmse = NULL
  )] %>% 
    
  # Bind pilots
  rbind(., pilot_preds) %>% 
  
  # Adjust for eta minus label
  .[, scen_summary := gsub(
    "eta1=-", 
    "eta1=min", 
    scen_summary
  )] %>%
  
  # Separate the scenario columns
  .[, (labs_scens) := data.table::tstrsplit(scen_summary, split = "-")] %>% 
  
  # Relevel methods
  .[, analy := factor(
    analy, levels = c("ref", "CCA", "ch1", "ch12", "ch12_int", "smcfcs")
  )] %>% 
  
  # Relevel other factors and properly names
  .[, ':=' (
    prop_miss = factor(
      prop_miss, 
      levels = c("prop_miss=0.1","prop_miss=0.5"),
      labels = c("10%", "50%")
    ),
    haz_shape = factor(
      haz_shape,
      levels = c("haz_shape=similar", "haz_shape=different"),
      labels = c("similar", "different")
    ),
    beta1 = factor(
      beta1,
      levels = c("beta1=0", "beta1=0.5", "beta1=1"),
      labels = c("0", "0.5", "1")
    ),
    eta1 = factor(
      eta1,
      levels = c("eta1=NA", "eta1=min1", "eta1=min2"),
      labels = c("None", "Weak", "Strong")
    ),
    miss_mech = factor(
      miss_mech,
      levels = c("miss_mech=MCAR", "miss_mech=MAR", 
                 "miss_mech=MAR_GEN", "miss_mech=MNAR"),
      labels = c("MCAR", "MAR", "MAR_GEN", "MNAR")
    ),
    X_level = factor(
      X_level,
      levels = c("X_level=continous", "X_level=binary"),
      labels = c("continuous", "binary")
    ),
    rho = factor(
      rho,
      levels = "rho=0.5",
      labels = "0.5"
    ),
    n = factor(
      n,
      levels = "n=2000",
      labels = "2000"
    )
  )] %>% 
  
  # Relevel prediction specific factors
  .[, ':=' (
    state = factor(
      state,
      levels = c("1", "2", "3"),
      labels = c("EFS", "REL", "NRM")
    ),
    times = factor(
      times, 
      levels = c(0.5, 5, 10),
      labels = c("6 months", "5 years", "10 years")
    )
  )]



# Export ready for plotting -----------------------------------------------


saveRDS(all_estims, file = "analysis/simulations/all_estims.rds")
saveRDS(all_preds, file = "analysis/simulations/all_preds.rds")


