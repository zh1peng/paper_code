#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=23:59:00
#SBATCH --job-name=per
#SBATCH --array=1-1000
#SBATCH --output=%x
cd /gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME_new_subid_unbalanced/PPI_L_VS
Rscript --vanilla LME_null_model.R