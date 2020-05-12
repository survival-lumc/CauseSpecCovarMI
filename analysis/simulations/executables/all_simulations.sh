#!/bin/bash
#SBATCH -J scenario_array
#SBATCH -n 1
#SBATCH --time=00:10:00 
#SBATCH --mem=1G
#SBATCH --partition=short
#SBATCH -o ./analysis/other/out_err_messages/scenario%a_job%A.out
#SBATCH --array=120-133 # These are the scenario numbers, 120-133

# run %2 so two running at any given times, each with 12 cores

# Set variables, change array to relevant scenario numbers
scenario=${SLURM_ARRAY_TASK_ID}

# Change scenarios to 14 (pilot) and replicate to 160
sbatch analysis/simulations/executables/one_simulation.sh $scenario 
	