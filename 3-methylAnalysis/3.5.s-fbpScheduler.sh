#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=32G
#SBATCH --job-name=3.5-featureBodyPlots
#SBATCH --output=/home/nia/project/nia/logs/%x

###################################################################
### SHELL SCHEDULER FOR PLOTTING METHYLATION ACROSS FEATURE BODIES
###################################################################

# load required modules
module restore rstudio_modules
# run script
Rscript 3.5-featureBodyPlots.R