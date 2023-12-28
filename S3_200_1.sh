#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M junming_shi@berkeley.edu
#$ -m base

export SIMU_NUM=3 SAMPLE_N=200 HALF1=TRUE HALF2=FALSE

R CMD BATCH --no-save run_simu/run_simu.R run_simu/simu_3_200_half1.Rout