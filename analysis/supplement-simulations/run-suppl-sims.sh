#!/bin/bash
#SBATCH -J suppl-sims
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:45:00 
#SBATCH --array=1-2

# Purge
module purge

# Load module
module load statistical/R/4.1.0/gcc.8.3.1

# Run imps - make different file for the other 
Rscript ./analysis/supplement-simulations/three-comprisks.R 
