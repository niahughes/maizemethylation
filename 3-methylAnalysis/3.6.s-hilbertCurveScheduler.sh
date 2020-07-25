#!/bin/bash
#SBATCH -t 9:00:00
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=3.6-hilbertCurve
#SBATCH --output=/home/nia/project/nia/logs/%x

###################################################################
### SHELL SCHEDULER FOR HILBERT CURVES OF DIFFERENTIAL METHYLATION
###################################################################

# load required modules
module restore rstudio_modules
# run script
Rscript 3.6-hilbertCurve.R