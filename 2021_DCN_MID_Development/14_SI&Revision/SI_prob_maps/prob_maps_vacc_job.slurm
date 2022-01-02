#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=prob_maps
#SBATCH --output=%x	
matlab -nodisplay -nosplash -nodesktop -r "run('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/prob_maps/run_prob_maps_job.m');"
