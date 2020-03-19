##********************************************##
## This is the file that gets called by Shark ##
##********************************************##

# Load scenarios 
scenarios <- readRDS("./inst/testdata/scenarios.rds") 

# Load all functions in compendium
devtools::load_all()

# Read-in those command line args
args <- commandArgs(TRUE)

# First bash argument is scenario, second is replicate
scen <- as.numeric(args[1])
repl <- as.numeric(args[2])

# Source in the one_simulation function
source("analysis/one_simulation.R")

# Run one replication of one scenario
one_simulation(scenario = scenarios[scenarios$scen_num == scen, ],
               rep_num = repl)


