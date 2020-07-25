#!/bin/bash
#SBATCH --time=8:00:00 
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=15G
#SBATCH --job-name=bam2nuc
#SBATCH --output=/scratch/nia/thesis/logs/%x

###################################################################
### BISMARK BAM2NUC: NUCLEOTIDE STATS ON ALIGNMENTS
###################################################################

# run bam2nuc on each alignment
bam2nuc --genome_folder /scratch/nia/thesis/genomes/A /scratch/nia/thesis/alignments/aMut/Amop1_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
bam2nuc --genome_folder /scratch/nia/thesis/genomes/A /scratch/nia/thesis/alignments/aWT/A_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
bam2nuc --genome_folder /scratch/nia/thesis/genomes/B /scratch/nia/thesis/alignments/bMut/Bmop1_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
bam2nuc --genome_folder /scratch/nia/thesis/genomes/B /scratch/nia/thesis/alignments/bWT/B_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
bam2nuc --genome_folder /scratch/nia/thesis/genomes/C /scratch/nia/thesis/alignments/cMut/Cmop1_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
bam2nuc --genome_folder /scratch/nia/thesis/genomes/C /scratch/nia/thesis/alignments/cWT/C_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
