#!/bin/bash
#SBATCH -J summarise_sims
#SBATCH -n 1
#SBATCH --time=00:30:00 
#SBATCH --mem=5G 
#SBATCH -o ~/job%A.out

# Purge
module purge

# Load module
module load statistical/R/4.0.2

# Run summary script
Rscript analysis/simulations/summarise-simulations.R
