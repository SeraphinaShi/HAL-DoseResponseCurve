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

    # Loop over combinations of HALF1 and HALF2
    for HALF1 in TRUE FALSE; do
      if [ "$HALF1" = "TRUE" ]; then
        HALF2=FALSE
        half_id=1
      else
        HALF2=TRUE
        half_id=2
      fi

      # Export the variables
      export SIMU_NUM SAMPLE_N HALF1 HALF2

      # Run the command
      R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_half${half_id}.Rout"

    done
  done
done