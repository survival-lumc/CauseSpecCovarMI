#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N onesim_MI
#$ -cwd
#$ -V
#$ -j Y
#$ -pe BWA 1
#$ -l h_vmem=10G

# SET VARIABLES
i=$1

R --vanilla --slave -f  one_simulation.R --args $i