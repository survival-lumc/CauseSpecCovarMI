#!/bin/bash

# Set variables
scenario=$i
replicate=$j

# Change scenarios to 14 (pilot) and replicate to 160
for scenario in {1..1}
do
	for replicate in {1..1}
	do
		sh analysis/one_simulation.sh $scenario $replicate
	done	
done
