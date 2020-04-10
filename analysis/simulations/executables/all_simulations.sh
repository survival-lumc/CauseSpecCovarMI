#!/bin/bash
#SBATCH -J scenario_array
#SBATCH -n 1
#SBATCH --time=00:05:00 
#SBATCH --mem=1G
#SBATCH --partition=short
#SBATCH -o ./analysis/other/out_err_messages/scenario%a_job%A.out
#SBATCH --array=15,16 # These are the scenario numbers

# run %2 so two running at any given times, each with 12 cores

# Set variables, change array to relevant scenario numbers
scenario=${SLURM_ARRAY_TASK_ID}

# Change scenarios to 14 (pilot) and replicate to 160
sbatch analysis/simulations/executables/one_simulation.sh $scenario 
	