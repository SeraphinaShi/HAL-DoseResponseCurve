#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu
#$ -m base

export SIMU_NUM=2 SAMPLE_N=500 HALF1=FALSE HALF2=FALSE GRID_EXTRA=TRUE

R CMD BATCH --no-save run_simu/run_simu.R run_simu/simu_2_500_gridextra.Rout