#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu


# Export the variables
export SIMU_NUM=1 SAMPLE_N=200 PART1=True PART2=False PART3=False
R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_part1.Rout"
