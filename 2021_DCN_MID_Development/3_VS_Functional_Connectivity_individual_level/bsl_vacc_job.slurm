#! /bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --job-name=bsl_ppi
#SBATCH --array=1-491
#SBATCH --output=%x_log_sub1k%a.out
cd $SLURM_SUBMIT_DIR	
matlab -r 'bsl_PPI_1st_level_vacc_function(${SLURM_ARRAY_TASK_ID}+1000)'
