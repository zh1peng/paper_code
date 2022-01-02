#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --job-name=bsl_merge
#SBATCH --output=%x
cd /gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/test_var/g2_merge_cons
python merge_con_files_by_subid.py