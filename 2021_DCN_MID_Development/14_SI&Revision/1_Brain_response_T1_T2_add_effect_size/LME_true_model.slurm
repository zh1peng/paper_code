#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=21
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --job-name=activiation_LME
#SBATCH --output=%x
cd /gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/activiation_LME_new_subid_unbalanced/add_effect_size
Rscript --vanilla LME_true_model.R