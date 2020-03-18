#!/bin/bash

# Set variables
scenario=$1 # First argument, the scenario
replicate=$2 # Second argument, the replicate

#echo "Scenario =" $i ", replication =" $j

Rscript analysis/run_simulations.R $scenario $replicate
 