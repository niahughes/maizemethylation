#!/bin/bash

###################################################################
### BISMARK ALIGNMENT DEDUPLICATION
###################################################################

# load required modules
module load bowtie2/2.3.4.1
module load samtools/1.9
	
# deduplicate aligned reads
sbatch --time=5:00:00 --mem=15G --job-name=aWTdedup --output=/scratch/nia/thesis/logs/%x \
	deduplicate_bismark -p --bam /scratch/nia/thesis/alignments/aWT/A_BC_R1_val_1_bismark_bt2_pe.bam

sbatch --time=5:00:00 --mem=15G --job-name=aMutdedup --output=/scratch/nia/thesis/logs/%x \
	deduplicate_bismark -p --bam /scratch/nia/thesis/alignments/aMut/Amop1_BC_R1_val_1_bismark_bt2_pe.bam

sbatch --time=5:00:00 --mem=15G --job-name=bWTdedup --output=/scratch/nia/thesis/logs/%x \
	deduplicate_bismark -p --bam /scratch/nia/thesis/alignments/bWT/B_BC_R1_val_1_bismark_bt2_pe.bam	

sbatch --time=5:00:00 --mem=15G --job-name=bMutdedup --output=/scratch/nia/thesis/logs/%x \
	deduplicate_bismark -p --bam /scratch/nia/thesis/alignments/bMut/Bmop1_BC_R1_val_1_bismark_bt2_pe.bam

sbatch --time=5:00:00 --mem=15G --job-name=cWTdedup --output=/scratch/nia/thesis/logs/%x \
	deduplicate_bismark -p --bam /scratch/nia/thesis/alignments/cWT/C_BC_R1_val_1_bismark_bt2_pe.bam

sbatch --time=5:00:00 --mem=15G --job-name=cMutdedup --output=/scratch/nia/thesis/logs/%x \
	deduplicate_bismark -p --bam /scratch/nia/thesis/alignments/cMut/Cmop1_BC_R1_val_1_bismark_bt2_pe.bam
