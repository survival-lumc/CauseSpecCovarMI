#!/bin/bash
#SBATCH -J repl_array
#SBATCH -n 1
#SBATCH --time=01:30:00 
#SBATCH --mem=4G 
#SBATCH -o ~/job%A_replicate%a.out
#SBATCH --array=1-400:5 

# Read in args
scenario=$1 # First argument, the scenario
replicate=${SLURM_ARRAY_TASK_ID} # Second argument, the replicate

# Stop if 160

# Purge
module purge

# Load module
module load statistical/R/3.6.2/gcc.8.3.1.lua

# The below is how to batch
i=$i

# Batch size and end of sequence - minus one since we only want a batch of x
batch_size=5
end=$(($replicate + $batch_size - 1)) # double brackets is called 'arithmetic expansion', need dollar sign before though


# Add an if statement for size of n_sim


# Begin sequence
seqo=$(seq $replicate $end)

for i in $seqo
do
  if [ $i -gt 400 ] # Stop if index is greater than (-gt) n_sim, 160 for n = 200 and 625 for n = 500
  then
      break
  fi
	Rscript analysis/simulations/executables/run_simulations.R $scenario $i
done