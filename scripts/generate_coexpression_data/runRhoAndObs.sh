#!/usr/bin/env bash
#SBATCH --job-name=rho_split_400
#SBATCH --partition=any_cpu
#SBATCH --array=1-11630

#for each ORF (one at a time) get its coexpression with all other ORFs and get the number of observations for each pair (ie how many samples where both ORFs are detetcted)
#repear for all 11630 ORFs
Rscript getRhoAndObs.R ${SLURM_ARRAY_TASK_ID} clr_raw5_sample400.RDS raw5_sample400

