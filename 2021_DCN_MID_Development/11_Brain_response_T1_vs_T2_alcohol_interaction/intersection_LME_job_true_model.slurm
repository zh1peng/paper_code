#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --job-name=activiation_LME
#SBATCH --output=%x
cd /gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/activiation_alc_LME_new_subid_balanced
Rscript --vanilla intersection_LME_true_model.R