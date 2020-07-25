#!/bin/bash
#SBATCH --time=5:00:00 
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --job-name=bedtoolsCoverage
#SBATCH --output=/home/nia/project/nia/logs/%x

###################################################################
### BEDTOOLS ALIGNMENT COVERAGE STATISTICS
###################################################################

# load required modules
module load bedtools

# get genome coverage stats in .txt form
bedtools genomecov -ibam /scratch/nia/thesis/alignments/aMut/Amop1_sorted.bam -max 75 > /scratch/nia/thesis/alignments/coverage/Amop1cov.txt
bedtools genomecov -ibam /scratch/nia/thesis/alignments/aWT/A_sorted.bam -max 75 > /scratch/nia/thesis/alignments/coverage/Acov.txt
bedtools genomecov -ibam /scratch/nia/thesis/alignments/bMut/Bmop1_sorted.bam -max 75 > /scratch/nia/thesis/alignments/coverage/Bmop1cov.txt
bedtools genomecov -ibam /scratch/nia/thesis/alignments/bWT/B_sorted.bam -max 75 > /scratch/nia/thesis/alignments/coverage/Bcov.txt
bedtools genomecov -ibam /scratch/nia/thesis/alignments/cMut/Cmop1_sorted.bam -max 75 > /scratch/nia/thesis/alignments/coverage/Cmop1cov.txt
bedtools genomecov -ibam /scratch/nia/thesis/alignments/cWT/C_sorted.bam -max 75 > /scratch/nia/thesis/alignments/coverage/Ccov.txt

# get genome coverage stats in .bedgraph form
bedtools genomecov -bga -ibam /home/nia/project/nia/alignments/aMut/Amop1_10.bam -max 10 > /home/nia/project/nia/alignments/coverage/Amop1cov.bedgraph
bedtools genomecov -bga -ibam /home/nia/project/nia/alignments/aWT/A_10.bam -max 10 > /home/nia/project/nia/alignments/coverage/Acov.bedgraph
bedtools genomecov -bga -ibam /home/nia/project/nia/alignments/bMut/Bmop1_10.bam -max 10 > /home/nia/project/nia/alignments/coverage/Bmop1cov.bedgraph
bedtools genomecov -bga -ibam /home/nia/project/nia/alignments/bWT/B_10.bam -max 10 > /home/nia/project/nia/alignments/coverage/Bcov.bedgraph
bedtools genomecov -bga -ibam /home/nia/project/nia/alignments/cMut/Cmop1_10.bam -max 10 > /home/nia/project/nia/alignments/coverage/Cmop1cov.bedgraph
bedtools genomecov -bga -ibam /home/nia/project/nia/alignments/cWT/C_10.bam -max 10 > /home/nia/project/nia/alignments/coverage/Ccov.bedgraph