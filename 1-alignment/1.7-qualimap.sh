#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=qualimap
#SBATCH --output=/home/nia/project/nia/logs/%x

###################################################################
### QUALIMAP ALIGNMENT QUALITY ASSESSMENT
###################################################################

# load required modules
module add mugqic/qualimap/2.2.1

# prepare cluster display for java commands
unset DISPLAY
echo $DISPLAY

# run qualimap individually on each subset, sorted, and indexed alignment file 
qualimap bamqc -bam /home/nia/project/nia/alignments/aMut/Amop1_10.bam --java-mem-size=5G -c -nt 8 -outdir /home/nia/project/nia/results/qualimap/aMut/
qualimap bamqc -bam /home/nia/project/nia/alignments/aWT/A_10.bam --java-mem-size=5G -c -nt 8 -outdir /home/nia/project/nia/results/qualimap/aWT/
qualimap bamqc -bam /home/nia/project/nia/alignments/bMut/Bmop1_10.bam --java-mem-size=5G -c -nt 8 -outformat HTML -outdir /home/nia/project/nia/results/qualimap/bMut/
qualimap bamqc -bam /home/nia/project/nia/alignments/bWT/B_10.bam --java-mem-size=5G -c -nt 8 -outdir /home/nia/project/nia/results/qualimap/bWT/
qualimap bamqc -bam /home/nia/project/nia/alignments/cMut/Cmop1_10.bam --java-mem-size=5G -c -nt 8 -outdir /home/nia/project/nia/results/qualimap/cMut/
qualimap bamqc -bam /home/nia/project/nia/alignments/cWT/C_10.bam --java-mem-size=5G -c -nt 8 -outdir /home/nia/project/nia/results/qualimap/cWT/

# run qualimap collectively on all alignment files
qualimap multi-bamqc -d 1.7.1-qualimapDir.txt -outdir /home/nia/project/nia/results/qualimap/all3/