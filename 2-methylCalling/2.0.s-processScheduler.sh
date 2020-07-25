#!/bin/bash
#SBATCH -t 10:00:00 
#SBATCH -N 6
#SBATCH -n 6
#SBATCH -c 1
#SBATCH --mem-per-cpu=10G
#SBATCH --job-name=2.0.1-processAll
#SBATCH --output=/home/nia/project/nia/logs/%x

###################################################################
### SHELL SCHEDULER FOR METHYLKIT PROCESSBISMARKALN SCRIPTS
### These scripts read in Bismark alignment files and 
### calculate methylation proportions for each cytosine, then
### output context-specific .txt files for use by 
### methylKit command methRead
###################################################################

# load required modules
module restore rstudio_modules

# run processBismarkAln on all files
srun -N 1 -n 1 Rscript 2.0-process_aWT.R &
srun -N 1 -n 1 Rscript 2.0-process_aMut.R &
srun -N 1 -n 1 Rscript 2.0-process_bWT.R &
srun -N 1 -n 1 Rscript 2.0-process_bMut.R &
srun -N 1 -n 1 Rscript 2.0-process_cWT.R &
srun -N 1 -n 1 Rscript 2.0-process_cMut.R &
wait