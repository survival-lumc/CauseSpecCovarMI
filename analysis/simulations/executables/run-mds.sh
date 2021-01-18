#!/bin/bash
#SBATCH -J mds-parall
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=03:00:00 
#SBATCH --mem=5G 


# Purge
module purge

# Load module
module load statistical/R/4.0.2

# Run imps
Rscript ./analysis/run-mds-imputations.R