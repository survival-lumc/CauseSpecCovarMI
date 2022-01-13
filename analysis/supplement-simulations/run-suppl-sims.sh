#!/bin/bash
#SBATCH -J suppl-sims
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:45:00 
#SBATCH --mem=5G 


# Purge
module purge

# Load module
module load statistical/R/4.1.0/gcc.8.3.1

# Run imps
Rscript ./analysis/supplement-simulations/three-comprisks.R