#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=2:59:00
#SBATCH --job-name=vertx

cd /gpfs1/home/z/c/zcao4/PCA_project/IMAGEN_FU2_vertex
Rscript --vanilla pca_vertex_cmd.R 25 lm