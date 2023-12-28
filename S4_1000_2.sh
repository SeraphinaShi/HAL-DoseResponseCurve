#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu
#$ -m base

export SIMU_NUM=4 SAMPLE_N=1000 HALF1=FALSE HALF2=TRUE

R CMD BATCH --no-save run_simu/run_simu.R run_simu/simu_4_1000_half2.Rout