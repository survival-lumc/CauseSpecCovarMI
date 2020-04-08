#!/bin/bash
#SBATCH -J repl_array
#SBATCH -n 1
#SBATCH --time=00:55:00 
#SBATCH --mem=4G 
#SBATCH --partition=short
#SBATCH -o ./analysis/other/out_err_messages/job%A_replicate%a.out
#SBATCH --array=1-160:3%36

# Value after % eg. %12 should be multiple of 4?
# %A is master task ID, and 
scenario=$1 # First argument, the scenario
replicate=${SLURM_ARRAY_TASK_ID} # Second argument, the replicate

#i=$i

# Purge
module purge

# Load module
module load statistical/R/3.6.2/gcc.8.3.1.lua

# Uncomment next line if batching happening:
#Rscript analysis/simulations/executables/run_simulations.R $scenario $replicate

# The below is how to batch - must be proceede by 1-10:3 (for batch size = 3)

i=$i

# Batch size and end of sequence - minus one since we only want a batch of x
batch_size=3
end=$(($replicate + $batch_size - 1)) # double brackets is called 'arithmetic expansion', need dollar sign before though

# Begin sequence
seqo=$(seq $replicate $end)

for i in $seqo
do
	Rscript analysis/simulations/executables/run_simulations.R $scenario $i
done