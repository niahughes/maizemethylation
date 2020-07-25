#!/bin/bash

###################################################################
### BISMARK ALIGNMENT TO CUSTOMIZED GENOMES
###################################################################

# load required modules
module load bowtie2/2.3.4.1
module load samtools/1.9
	
# aLign to respective genomes
sbatch --time=72:00:00 --nodes=1 --ntasks-per-node=32 --mem=0 --job-name=aWTalignment --output=/scratch/nia/thesis/logs/%x \
	bismark --bowtie2 --multicore 8 -p 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "A_BC" \
	--genome /scratch/nia/thesis/genomes/A/ \
	--temp_dir /scratch/nia/thesis/temp/aWT/ \
	-o /scratch/nia/thesis/alignments/aWT/ \
	-1 /scratch/nia/thesis/reads/A_BC_R1_val_1.fq.gz \
	-2 /scratch/nia/thesis/reads/A_BC_R2_val_2.fq.gz

sbatch --time=72:00:00 --nodes=1 --ntasks-per-node=32 --mem=0 --job-name=aMutalignment --output=/scratch/nia/thesis/logs/%x \
	bismark --bowtie2 --multicore 8 -p 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "Amop1_BC" \
	--genome /scratch/nia/thesis/genomes/A/ \
	--temp_dir /scratch/nia/thesis/temp/aMut/ \
	-o /scratch/nia/thesis/alignments/aMut/ \
	-1 /scratch/nia/thesis/reads/Amop1_BC_R1_val_1.fq.gz \
	-2 /scratch/nia/thesis/reads/Amop1_BC_R2_val_2.fq.gz

sbatch --time=72:00:00 --nodes=1 --ntasks-per-node=32 --mem=0 --job-name=bWTalignment --output=/scratch/nia/thesis/logs/%x \
	bismark --bowtie2 --multicore 8 -p 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "B_BC" \
	--genome /scratch/nia/thesis/genomes/B/ \
	--temp_dir /scratch/nia/thesis/temp/bWT/ \
	-o /scratch/nia/thesis/alignments/bWT/ \
	-1 /scratch/nia/thesis/reads/B_BC_R1_val_1.fq.gz \
	-2 /scratch/nia/thesis/reads/B_BC_R2_val_2.fq.gz

sbatch --time=72:00:00 --nodes=1 --ntasks-per-node=32 --mem=0 --job-name=bMutalignment --output=/scratch/nia/thesis/logs/%x \
	bismark --bowtie2 --multicore 8 -p 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "Bmop1_BC" \
	--genome /scratch/nia/thesis/genomes/B/ \
	--temp_dir /scratch/nia/thesis/temp/bMut/ \
	-o /scratch/nia/thesis/alignments/bMut/ \
	-1 /scratch/nia/thesis/reads/Bmop1_BC_R1_val_1.fq.gz \
	-2 /scratch/nia/thesis/reads/Bmop1_BC_R2_val_2.fq.gz

sbatch --time=72:00:00 --nodes=1 --ntasks-per-node=32 --mem=0 --job-name=cWTalignment --output=/scratch/nia/thesis/logs/%x \
	bismark --bowtie2 --multicore 8 -p 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "C_BC" \
	--genome /scratch/nia/thesis/genomes/C/ \
	--temp_dir /scratch/nia/thesis/temp/cWT/ \
	-o /scratch/nia/thesis/alignments/cWT/ \
	-1 /scratch/nia/thesis/reads/C_BC_R1_val_1.fq.gz \
	-2 /scratch/nia/thesis/reads/C_BC_R2_val_2.fq.gz

sbatch --time=72:00:00 --nodes=1 --ntasks-per-node=32 --mem=0 --job-name=cMutalignment --output=/scratch/nia/thesis/logs/%x \
	bismark --bowtie2 --multicore 8 -p 4 \
	-D 20 -R 3 -N 1 \
	--score_min L,0,-0.6 \
	--rg_tag --rg_id "1" --rg_sample "Cmop1_BC" \
	--genome /scratch/nia/thesis/genomes/C/ \
	--temp_dir /scratch/nia/thesis/temp/cMut/ \
	-o /scratch/nia/thesis/alignments/cMut/ \
	-1 /scratch/nia/thesis/reads/Cmop1_BC_R1_val_1.fq.gz \
	-2 /scratch/nia/thesis/reads/Cmop1_BC_R2_val_2.fq.gz