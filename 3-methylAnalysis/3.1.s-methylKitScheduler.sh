#!/bin/bash
#SBATCH --time=12:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-user=nia@uoguelph.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name=3.1.1-mKbedgraphs
#SBATCH --output=/home/nia/project/nia/logs/%x

###################################################################
### SHELL SCHEDULER FOR METHYLKIT BEDGRAPH GENERATION
###################################################################

# load required modules
module restore rstudio_modules
# run script
Rscript 3.1.1-mKbedgraphs.R

# sort bedgraphs by position
cd /home/nia/project/nia/results/bedgraphs/mKcoverage/
FILES=*
for f in $FILES;
do sort -k1,1 -k2,2n -n $f -o sorted.$f;
done

cd /home/nia/project/nia/results/bedgraphs/mKpercMeth/
FILES=*
for f in $FILES;
	do sort -k1,1 -k2,2n -n $f -o sorted.$f;
done

cd /home/nia/project/nia/results/bedgraphs/mKmethDiff/
FILES=*
for f in $FILES;
do sort -k1,1 -k2,2n -n $f -o sorted.$f;
done

# average perc.meth bedgraphs by mutation status
module load bedtools
cd /home/nia/project/nia/results/bedgraphs/mKpercMeth

bedtools unionbedg -filler NA -i sorted.aMutmKCpGpercMeth3x.bedgraph sorted.bMutmKCpGpercMeth3x.bedgraph sorted.cMutmKCpGpercMeth3x.bedgraph | awk '{sum=cnt=0; for (i=4;i<=NF;i++) if ($i != "NA") { sum+=$i; cnt++ } print $1"\t"$2"\t"$3"\t", (cnt ? sum/cnt : "NA") }' > poolCpGMut.bedgraph
bedtools unionbedg -filler NA -i sorted.aWTmKCpGpercMeth3x.bedgraph sorted.bWTmKCpGpercMeth3x.bedgraph sorted.cWTmKCpGpercMeth3x.bedgraph | awk '{sum=cnt=0; for (i=4;i<=NF;i++) if ($i != "NA") { sum+=$i; cnt++ } print $1"\t"$2"\t"$3"\t", (cnt ? sum/cnt : "NA") }' > poolCpGWT.bedgraph

bedtools unionbedg -filler NA -i sorted.aMutmKCHGpercMeth3x.bedgraph sorted.bMutmKCHGpercMeth3x.bedgraph sorted.cMutmKCHGpercMeth3x.bedgraph | awk '{sum=cnt=0; for (i=4;i<=NF;i++) if ($i != "NA") { sum+=$i; cnt++ } print $1"\t"$2"\t"$3"\t", (cnt ? sum/cnt : "NA") }' > poolCHGMut.bedgraph
bedtools unionbedg -filler NA -i sorted.aWTmKCHGpercMeth3x.bedgraph sorted.bWTmKCHGpercMeth3x.bedgraph sorted.cWTmKCHGpercMeth3x.bedgraph | awk '{sum=cnt=0; for (i=4;i<=NF;i++) if ($i != "NA") { sum+=$i; cnt++ } print $1"\t"$2"\t"$3"\t", (cnt ? sum/cnt : "NA") }' > poolCHGWT.bedgraph

bedtools unionbedg -filler NA -i sorted.aMutmKCHHpercMeth3x.bedgraph sorted.bMutmKCHHpercMeth3x.bedgraph sorted.cMutmKCHHpercMeth3x.bedgraph | awk '{sum=cnt=0; for (i=4;i<=NF;i++) if ($i != "NA") { sum+=$i; cnt++ } print $1"\t"$2"\t"$3"\t", (cnt ? sum/cnt : "NA") }' > poolCHHMut.bedgraph
bedtools unionbedg -filler NA -i sorted.aWTmKCHHpercMeth3x.bedgraph sorted.bWTmKCHHpercMeth3x.bedgraph sorted.cWTmKCHHpercMeth3x.bedgraph | awk '{sum=cnt=0; for (i=4;i<=NF;i++) if ($i != "NA") { sum+=$i; cnt++ } print $1"\t"$2"\t"$3"\t", (cnt ? sum/cnt : "NA") }' > poolCHHWT.bedgraph
