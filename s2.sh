#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu
#$ -m base
# Loop over SIMU_NUM
#!/bin/bash

# Start scripts in parallel
export SIMU_NUM=2 SAMPLE_N=200 PART1=True PART2=False PART3=False; R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_part1.Rout" &
export SIMU_NUM=2 SAMPLE_N=500 PART1=True PART2=False PART3=False; R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_part1.Rout" &
export SIMU_NUM=2 SAMPLE_N=1000 PART1=True PART2=False PART3=False; R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_part1.Rout" &

export SIMU_NUM=2 SAMPLE_N=200 PART1=False PART2=True PART3=False; R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_part2.Rout" &
export SIMU_NUM=2 SAMPLE_N=500 PART1=False PART2=True PART3=False; R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_part2.Rout" &
export SIMU_NUM=2 SAMPLE_N=1000 PART1=False PART2=True PART3=False; R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_part2.Rout" &

export SIMU_NUM=2 SAMPLE_N=200 PART1=False PART2=False PART3=True; R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_part3.Rout" &
export SIMU_NUM=2 SAMPLE_N=500 PART1=False PART2=False PART3=True; R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_part3.Rout" &
export SIMU_NUM=2 SAMPLE_N=1000 PART1=False PART2=False PART3=True; R CMD BATCH --no-save run_simu/run_simu.R "run_simu/simu_${SIMU_NUM}_${SAMPLE_N}_part3.Rout" &

# Wait for all background jobs to finish
wait