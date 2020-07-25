#!/bin/bash
#SBATCH --time=5:00:00 
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --job-name=fastqc
#SBATCH --output=/home/nia/project/nia/logs/%x

###################################################################
### QUALITY ANALYSIS OF READS
###################################################################

# load fastqc
module load fastqc/0.11.8 

# run fastqc on all files in reads folder
fastqc -o /home/nia/project/nia/results/fastqc/ /home/nia/project/nia/reads/* 