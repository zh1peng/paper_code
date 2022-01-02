#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=23:59:00
#SBATCH --job-name=per_R
#SBATCH --array=1-1000
#SBATCH --output=%x
cd /gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME_alc_new_subid_balanced/PPI_R_VS
Rscript --vanilla intersection_LME_null_model.R