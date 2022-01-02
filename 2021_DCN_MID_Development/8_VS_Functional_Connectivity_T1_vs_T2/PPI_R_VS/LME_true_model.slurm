#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=21
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --job-name=PPI_R_LME
#SBATCH --output=%x
cd /gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME_new_subid_unbalanced/PPI_R_VS
Rscript --vanilla LME_true_model.R