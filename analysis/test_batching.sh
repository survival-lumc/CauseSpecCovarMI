#!/bin/bash

# Set variables
scenario=$i

#index=$(seq 1 3 160) # This is the job array
batch_size=3
begin=7 # numbers do not need bracket, this is the #slurm_array_index
end=$(($begin + $batch_size -1)) # double brackets is called 'arithmetic expansion', need dollar sign before though
# Minus one since we only want a batch of 3
seqo=$(seq $begin $end)


for scenario in $seqo
do
	echo $scenario
done
