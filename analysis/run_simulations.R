##********************************************##
## This is the file that gets called by Shark ##
##********************************************##

# Load scenarios 
scenarios <- readRDS("inst/testdata/scenarios.rds") 

# Load all functions in compendium
devtools::load_all()

# Read-in those command line args
args <- commandArgs(TRUE)

# First is scenario, second is replicate
scen <- as.numeric(args[1])
repl <- as.numeric(args[2])

# Source in the one_simulation function
source("analysis/one_simulation.R")


one_simulation(scenario = scenarios[scenarios$scen_num == scen, ],
               rep_num = repl)


#View(test_run$estims %>%  dplyr::mutate_if(is.numeric, ~ round(., 3)))
#View(test_run$preds %>%  dplyr::mutate_if(is.numeric, ~ round(., 3))

# do the timings!