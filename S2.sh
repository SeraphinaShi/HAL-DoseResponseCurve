#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu

for SAMPLE_N in 200 500 1000; do

    # Export the variables
    export SIMU_NUM=2 SAMPLE_N
    
    R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_first_smallKnots.Rout"

done