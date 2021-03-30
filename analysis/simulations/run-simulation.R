##************************************##
## File called to run one replication ##
##************************************##

# Load scenarios 
scenarios_df <- CauseSpecCovarMI::scenarios

# Load all functions in compendium
library(CauseSpecCovarMI)

# Read-in those command line args
args <- commandArgs(TRUE)

# First bash argument is scenario, second is replicate
scen <- as.numeric(args[1])
repl <- as.numeric(args[2])

# Run one replication of one scenario
CauseSpecCovarMI::one_simulation(
  scenario = scenarios_df[scenarios_df$scen_num == scen, ],
  rep_num = repl
)

