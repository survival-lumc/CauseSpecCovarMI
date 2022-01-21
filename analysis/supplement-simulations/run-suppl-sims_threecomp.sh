#!/bin/bash
#SBATCH -J suppl-sims-threecomp
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=15:00:00 
#SBATCH --array=3-16 

# Purge
module purge

# Load module
module load statistical/R/4.1.0/gcc.8.3.1

# Run imps - make different file for the other 
Rscript ./analysis/supplement-simulations/three-comprisks.R 
