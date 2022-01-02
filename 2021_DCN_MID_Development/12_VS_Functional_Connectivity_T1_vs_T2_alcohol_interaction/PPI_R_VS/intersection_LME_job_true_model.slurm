#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --job-name=PPI_alc_R
#SBATCH --output=%x
cd /gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME_alc_new_subid_balanced/PPI_R_VS
Rscript --vanilla intersection_LME_true_model.R