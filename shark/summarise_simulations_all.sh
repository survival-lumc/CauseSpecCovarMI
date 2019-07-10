#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N summarise_sims
#$ -l h_vmem=5G
#$ -cwd
#$ -V
#$ -j Y
#$ -m be

module add R R/3.5.3

Rscript /exports/molepi/users/efbonneville/mi_comprisks/results_ISCB/summarise_simulations_all.R