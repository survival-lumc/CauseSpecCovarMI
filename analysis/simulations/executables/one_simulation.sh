#!/bin/bash
#SBATCH -J sim_causespec
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=01:00:00 
#SBATCH --mem=4G
#SBATCH -o ./out_err_messages/job_%J.out

# Set variables
scenario=$1 # First argument, the scenario
replicate=$2 # Second argument, the replicate

# Load module
module load statistical/R/3.6.2/gcc.8.3.1.lua

Rscript analysis/run_simulations.R $scenario $replicate
 