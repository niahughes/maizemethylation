#!/bin/bash

###################################################################
### BISMARK INDEXING OF CUSTOMIZED GENOMES
###################################################################

# load bowtie2
module load bowtie2/2.3.4.1

# index genomes
sbatch --time=8:00:00 --mem-per-cpu=4000 --ntasks=10 --job-name=aIndexBismark --output=/scratch/nia/thesis/logs/%x \
	bismark_genome_preparation --bowtie2 /scratch/nia/thesis/genomes/A/

sbatch --time=8:00:00 --mem-per-cpu=4000 --ntasks=10 --job-name=bIndexBismark --output=/scratch/nia/thesis/logs/%x \
	bismark_genome_preparation --bowtie2 /scratch/nia/thesis/genomes/B/

sbatch --time=8:00:00 --mem-per-cpu=4000 --ntasks=10 --job-name=cIndexBismark --output=/scratch/nia/thesis/logs/%x \
	bismark_genome_preparation --bowtie2 /scratch/nia/thesis/genomes/C/

# index chloroplast genome
sbatch --time=1:00:00 --mem-per-cpu=1G --ntasks=1 --job-name=indexPt --output=/home/nia/project/nia/logs/%x \
	bismark_genome_preparation --bowtie2 /scratch/nia/thesis/genomes/Pt/