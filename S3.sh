#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu

SIMU_NUM = 3

for SAMPLE_N in 200 500 1000; do

    HALF1 = FALSE
    HALF2 = FALSE
    GRID_EXTRA = FALSE
    FIRST_SMALLKNOTS = TRUE

    # Export the variables
    export SIMU_NUM SAMPLE_N HALF1 HALF2 GRID_EXTRA FIRST_SMALLKNOTS
    
    R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_first_smallKnots.Rout"

done