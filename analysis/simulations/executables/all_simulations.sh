#!/bin/bash
#SBATCH -J sim_causespec_head
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH --time=00:15:00 
#SBATCH -o ./out_err_messages/job_%J.out

# Set variables
scenario=$i
replicate=$j

# Change scenarios to 14 (pilot) and replicate to 160
for scenario in {3..3}
do
	for replicate in {1..1}
	do
		sbatch analysis/one_simulation.sh $scenario $replicate
	done	
done
