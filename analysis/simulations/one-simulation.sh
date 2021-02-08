#!/bin/bash
#SBATCH -J repl_array
#SBATCH -n 1
#SBATCH --time=01:30:00 
#SBATCH --mem=4G 
#SBATCH -o ~/job%A_replicate%a.out
#SBATCH --array=1-160:3 # Batch each job by 3

# Read in args
scenario=$1 # First argument, the scenario
replicate=${SLURM_ARRAY_TASK_ID} # Second argument, the replicate

# Purge
module purge

# Load module
module load statistical/R/4.0.2

# Batching
batch_size=3
end=$(($replicate + $batch_size - 1)) 
seqo=$(seq $replicate $end)

# Run each replication of batch (e.g. first batch is reps 1, 2 and 3
for rep in $seqo
do
  if [ $rep -gt 160 ] 
  then
	break
  fi
	Rscript analysis/simulations/run-simulation.R $scenario $rep
done