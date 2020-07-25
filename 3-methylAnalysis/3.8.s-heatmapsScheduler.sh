#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=32G
#SBATCH --job-name=3.8-heatmaps
#SBATCH --output=/home/nia/project/nia/logs/%x

###################################################################
### SHELL SCHEDULER FOR FEATURE-SPECIFIC METHYLATION CHANGE 
### HEATMAPS
###################################################################

# load required modules
module restore rstudio_modules
# run script
Rscript 3.8-heatmaps.R