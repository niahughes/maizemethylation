#!/bin/bash
#SBATCH --time=5:00:00 
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=subsetBAM
#SBATCH --output=/home/nia/project/nia/logs/%x

###################################################################
### SUBSET ALIGNMENTS TO READS MAPPING TO CHROMOSOMES 1-10
###################################################################

# subset alignments
samtools view /home/nia/project/nia/alignments/aMut/Emop1_sorted.bam {1..10} -b -o /home/nia/project/nia/alignments/aMut/Amop1_10.bam
samtools view /home/nia/project/nia/alignments/aWT/E_sorted.bam {1..10} -b -o /home/nia/project/nia/alignments/aWT/A_10.bam
samtools view /home/nia/project/nia/alignments/bMut/Gmop1_sorted.bam {1..10} -b -o /home/nia/project/nia/alignments/bMut/Bmop1_10.bam
samtools view /home/nia/project/nia/alignments/bWT/G_sorted.bam {1..10} -b -o /home/nia/project/nia/alignments/bWT/B_10.bam
samtools view /home/nia/project/nia/alignments/cMut/Hmop1_sorted.bam {1..10} -b -o /home/nia/project/nia/alignments/cMut/Cmop1_10.bam
samtools view /home/nia/project/nia/alignments/cWT/H_sorted.bam {1..10} -b -o /home/nia/project/nia/alignments/cWT/C_10.bam

# index new subset BAM files
samtools index /home/nia/project/nia/alignments/aMut/Amop1_10.bam
samtools index /home/nia/project/nia/alignments/aWT/A_10.bam
samtools index /home/nia/project/nia/alignments/bMut/Bmop1_10.bam
samtools index /home/nia/project/nia/alignments/bWT/B_10.bam
samtools index /home/nia/project/nia/alignments/cMut/Cmop1_10.bam
samtools index /home/nia/project/nia/alignments/cWT/C_10.bam