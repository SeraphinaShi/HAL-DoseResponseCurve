#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu
#$ -m base
# Loop over SIMU_NUM
for SIMU_NUM in 1 2 3 4; do

  # Loop over SAMPLE_N
  for SAMPLE_N in 200 500 1000; do

    HALF1 = TRUE
    HALF2 = TRUE
    GRID_EXTRA = TRUE
    FIRST_SMALLKNOTS = TRUE

    # Export the variables
    export SIMU_NUM SAMPLE_N HALF1 HALF2 GRID_EXTRA FIRST_SMALLKNOTS
    
    R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}.Rout"

  done
done