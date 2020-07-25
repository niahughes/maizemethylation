#!/bin/bash
#SBATCH --time=24:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=24G
#SBATCH --job-name=2.2-databases
#SBATCH --output=/home/nia/project/nia/logs/%x

###MethylKit analysis: database creation scheduler
##3/14/2020

###################################################################
### SHELL SCHEDULER FOR METHYLKIT METHYLATION CALLING
### This script reads filtered tabix files of methylation 
### proportions, tiles into 100bp windows, removes introgression-
### affected region from chromosome 2, unites methylation
### calls, calculates differential methylation, and calls 
### significantly different tiles for each context
###################################################################

# load required modules
module restore rstudio_modules

# run script
Rscript 2.2-databases.R
