#!/bin/bash

###################################################################
### BISMARK EXTRACTION OF METHYLATION CALLS
###################################################################

# load required modules
module swap gcc intel
module load bowtie2/2.3.4.1
module load samtools/1.9

# run bismark methylation extractor on each deduplicated alignment
sbatch --time=5:00:00 --mem-per-cpu=5G --ntasks=15 --job-name=aWTextract --output=/scratch/nia/thesis/logs/%x \
	bismark_methylation_extractor  -p --comprehensive --multicore 5 -o /scratch/nia/thesis/alignments/aWT/extractor \
	/scratch/nia/thesis/alignments/aWT/A_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam

sbatch --time=5:00:00 --mem-per-cpu=5G --ntasks=15 --job-name=aMutextract --output=/scratch/nia/thesis/logs/%x \
	bismark_methylation_extractor  -p --comprehensive --multicore 5 -o /scratch/nia/thesis/alignments/aMut/extractor \
	/scratch/nia/thesis/alignments/aMut/Amop1_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam

sbatch --time=5:00:00 --mem-per-cpu=5G --ntasks=15 --job-name=bWTextract --output=/scratch/nia/thesis/logs/%x \
	bismark_methylation_extractor  -p --comprehensive --multicore 5 -o /scratch/nia/thesis/alignments/bWT/extractor \
	/scratch/nia/thesis/alignments/bWT/B_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam

sbatch --time=5:00:00 --mem-per-cpu=5G --ntasks=15 --job-name=bMutextract --output=/scratch/nia/thesis/logs/%x \
	bismark_methylation_extractor  -p --comprehensive --multicore 5 -o /scratch/nia/thesis/alignments/bMut/extractor \
	/scratch/nia/thesis/alignments/bMut/Bmop1_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam
	
sbatch --time=5:00:00 --mem-per-cpu=5G --ntasks=15 --job-name=cWTextract --output=/scratch/nia/thesis/logs/%x \
	bismark_methylation_extractor  -p --comprehensive --multicore 5 -o /scratch/nia/thesis/alignments/cWT/extractor \
	/scratch/nia/thesis/alignments/cWT/C_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam

sbatch --time=5:00:00 --mem-per-cpu=5G --ntasks=15 --job-name=cMutextract --output=/scratch/nia/thesis/logs/%x \
	bismark_methylation_extractor  -p --comprehensive --multicore 5 -o /scratch/nia/thesis/alignments/cMut/extractor \
	/scratch/nia/thesis/alignments/cMut/Cmop1_BC_R1_val_1_bismark_bt2_pe.deduplicated.bam