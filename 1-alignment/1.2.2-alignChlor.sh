#!/bin/bash

###################################################################
### BISMARK ALIGNMENT TO CHLOROPLAST GENOME TO DETERMINE FDR
###################################################################

# load required modules
module load bowtie2/2.3.4.1
module load samtools/1.9

# align reads to chloroplast genome
sbatch --time=8:00:00 --ntasks=4 --mem-per-cpu=5G --job-name=aWTchloroplast --output=/home/nia/project/nia/logs/FDR/%x \
	bismark --bowtie2  \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "A_BC" \
	--genome /home/nia/project/nia/genomes/Pt/ \
	--temp_dir /scratch/nia/thesis/temp/aWT/ \
	-o /home/nia/project/nia/alignments/FDR/ \
	-1 /home/nia/project/nia/reads/A_BC_R1_val_1.fq.gz \
	-2 /home/nia/project/nia/reads/A_BC_R2_val_2.fq.gz
	
sbatch --time=8:00:00 --ntasks=12 --mem-per-cpu=10G --job-name=aMutchloroplast --output=/home/nia/project/nia/logs/FDR/%x \
	bismark --bowtie2 --multicore 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "Amop1_BC" \
	--genome /home/nia/project/nia/genomes/Pt/ \
	--temp_dir /scratch/nia/thesis/temp/aMut/ \
	-o /home/nia/project/nia/alignments/FDR/ \
	-1 /home/nia/project/nia/reads/Amop1_BC_R1_val_1.fq.gz \
	-2 /home/nia/project/nia/reads/Amop1_BC_R2_val_2.fq.gz
	
sbatch --time=8:00:00 --ntasks=12 --mem-per-cpu=10G --job-name=bWTchloroplast --output=/home/nia/project/nia/logs/FDR/%x \
	bismark --bowtie2 --multicore 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "B_BC" \
	--genome /home/nia/project/nia/genomes/Pt/ \
	--temp_dir /scratch/nia/thesis/temp/bWT/ \
	-o /home/nia/project/nia/alignments/FDR/ \
	-1 /home/nia/project/nia/reads/B_BC_R1_val_1.fq.gz \
	-2 /home/nia/project/nia/reads/B_BC_R2_val_2.fq.gz
	
sbatch --time=8:00:00 --ntasks=12 --mem-per-cpu=10G --job-name=bMutchloroplast --output=/home/nia/project/nia/logs/FDR/%x \
	bismark --bowtie2 --multicore 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "Bmop1_BC" \
	--genome /home/nia/project/nia/genomes/Pt/ \
	--temp_dir /scratch/nia/thesis/temp/bMut/ \
	-o /home/nia/project/nia/alignments/FDR/ \
	-1 /home/nia/project/nia/reads/Bmop1_BC_R1_val_1.fq.gz \
	-2 /home/nia/project/nia/reads/Bmop1_BC_R2_val_2.fq.gz
	
sbatch --time=8:00:00 --ntasks=12 --mem-per-cpu=10G --job-name=cWTchloroplast --output=/home/nia/project/nia/logs/FDR/%x \
	bismark --bowtie2 --multicore 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "C_BC" \
	--genome /home/nia/project/nia/genomes/Pt/ \
	--temp_dir /scratch/nia/thesis/temp/cWT/ \
	-o /home/nia/project/nia/alignments/FDR/ \
	-1 /home/nia/project/nia/reads/C_BC_R1_val_1.fq.gz \
	-2 /home/nia/project/nia/reads/C_BC_R2_val_2.fq.gz
	
sbatch --time=8:00:00 --ntasks=12 --mem-per-cpu=10G --job-name=cMutchloroplast --output=/home/nia/project/nia/logs/FDR/%x \
	bismark --bowtie2 --multicore 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "Cmop1_BC" \
	--genome /home/nia/project/nia/genomes/Pt/ \
	--temp_dir /scratch/nia/thesis/temp/cMut/ \
	-o /home/nia/project/nia/alignments/FDR/ \
	-1 /home/nia/project/nia/reads/Cmop1_BC_R1_val_1.fq.gz \
	-2 /home/nia/project/nia/reads/Cmop1_BC_R2_val_2.fq.gz