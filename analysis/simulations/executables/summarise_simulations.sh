#!/bin/bash
#SBATCH -J summarise_sims
#SBATCH -n 1
#SBATCH --time=00:15:00 
#SBATCH --mem=10G
#SBATCH --partition=short
#SBATCH -o ./analysis/other/out_err_messages/job_%J.out


# Load module
module load statistical/R/3.6.2/gcc.8.3.1.lua

Rscript analysis/simulations/executables/summarise_simulations.R 
 