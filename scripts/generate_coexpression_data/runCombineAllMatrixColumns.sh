#!/usr/bin/env bash
#SBATCH --job-name=combine_files
#SBATCH --partition=dept_cpu

Rscript combineAllMatrixColumns.R ./rho rho_raw5_sample400.RDS
Rscript combineAllMatrixColumns.R ./num_obs numobs_raw5_sample400.RDS
