#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=32G
#SBATCH --job-name=2.1-methReadandFilter
#SBATCH --output=/home/nia/project/nia/logs/%x

###################################################################
### SHELL SCHEDULER FOR METHYLKIT METHREAD AND FILTERING
### This script reads in .txt files of methylation proportions,
### creates tabix databases of methylKit objects, and filters
### based on coverage (minimum 3 reads at each site, removal of 
### sites within top 0.1% of reads)
###################################################################

# load required modules
module restore rstudio_modules

# run script
Rscript 2.1-methReadandFilter.R