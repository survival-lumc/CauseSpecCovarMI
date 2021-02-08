#!/bin/bash
#SBATCH -J scenario_array
#SBATCH -n 1
#SBATCH --time=00:15:00 
#SBATCH --mem=1G
#SBATCH --partition=short # Change to appropriate partition
#SBATCH -o ~/scenario%a_job%A.out
#SBATCH --array=1-119 # These are the scenario numbers

# Set scenario variable
scenario=${SLURM_ARRAY_TASK_ID}

# Call all replication for single scenario
sbatch ./analysis/simulations/one-simulation.sh $scenario 
	