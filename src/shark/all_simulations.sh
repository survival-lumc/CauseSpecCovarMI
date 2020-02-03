#!/bin/bash
#$ -S /bin/bash
#$ -N sims_scenarios
#$ -q all.q
#$ -cwd
#$ -V
#$ -j Y

module add R R/3.5.3

# Change 1000 to however many datasets need for all scenarios
for i in {1..1000}
    do
	    qsub one_simulation.sh $i
    done