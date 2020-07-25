#!/bin/bash
#SBATCH --time=5:00:00 
#SBATCH --ntasks=5
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=sortIndexAlignments
#SBATCH --output=/scratch/nia/thesis/logs/%x

###################################################################
### SORT AND INDEX BISMARK OUTPUT FOR METHYLKIT IMPORT
###################################################################

# load required modules
module load samtools/1.9

# sort Bismark alignment files
samtools sort -o /scratch/nia/thesis/alignments/aWT/E_sorted.bam -@ 5 /scratch/nia/thesis/alignments/aWT/A_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
samtools sort -o /scratch/nia/thesis/alignments/aMut/Emop1_sorted.bam -@ 5 /scratch/nia/thesis/alignments/aMut/Amop1_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
samtools sort -o /scratch/nia/thesis/alignments/bWT/G_sorted.bam -@ 5 /scratch/nia/thesis/alignments/bWT/B_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
samtools sort -o /scratch/nia/thesis/alignments/bMut/Gmop1_sorted.bam -@ 5 /scratch/nia/thesis/alignments/bMut/Bmop1_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
samtools sort -o /scratch/nia/thesis/alignments/cWT/H_sorted.bam -@ 5 /scratch/nia/thesis/alignments/cWT/C_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
samtools sort -o /scratch/nia/thesis/alignments/cMut/Hmop1_sorted.bam -@ 5 /scratch/nia/thesis/alignments/cMut/Cmop1_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam

# index sorted Bismark alignments
samtools index /scratch/nia/thesis/alignments/aWT/A_sorted.bam
samtools index /scratch/nia/thesis/alignments/aMut/Amop1_sorted.bam
samtools index /scratch/nia/thesis/alignments/bWT/B_sorted.bam
samtools index /scratch/nia/thesis/alignments/bMut/Bmop1_sorted.bam
samtools index /scratch/nia/thesis/alignments/cWT/C_sorted.bam
samtools index /scratch/nia/thesis/alignments/cMut/Cmop1_sorted.bam