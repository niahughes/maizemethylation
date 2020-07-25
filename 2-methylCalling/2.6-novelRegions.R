###################################################################
### LONG DMRs NOT MAPPING TO FEATURES
### R script to find DMRs >= 400bp not mapping to any known gene
### or transposable element
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(ComplexHeatmap)

# get annotation files
geneObj <- readBed("/home/nia/project/nia/annotation/geneAnnotationb73.bed", remove.unusual = TRUE)
up <- flank(geneObj, width=2000, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE) # get upstream flank of feature region
down <- flank(geneObj, width=2000, start=FALSE, both=FALSE, use.names=TRUE, ignore.strand=FALSE) # get downstream flank of feature region
te.file <- keepStandardChromosomes(import.gff("/home/nia/project/nia/annotation/B73.structuralTEv2.fulllength.2018-09-19.gff3", format="gff3"), pruning.mode="coarse")
legacy.file <- keepStandardChromosomes(import.gff("/home/nia/project/nia/annotation/legacy.gff3", format="gff3"), pruning.mode="coarse")
# get reference genome sequence
refGenome <- readDNAStringSet("/home/nia/project/nia/genomes/Zea_mays.B73_RefGen_v4.dna.toplevel.fa")
# list DM tiles for each comparison
inFiles <- list(fullCpG=CpG.3x.100.oChi.50p.all, aCpG=aCpG.3x.100.50p.all, bCpG=bCpG.3x.100.50p.all, cCpG=cCpG.3x.100.50p.all, # note no bCHH - no tiles match results
                fullCHG=CHG.3x.100.oChi.50p.all, aCHG=aCHG.3x.100.50p.all, bCHG=bCHG.3x.100.50p.all, hCHG=hCHG.3x.100.50p.all,
                fullCHH=CHH.3x.100.oChi.50p.all, aCHH=aCHH.3x.100.50p.all, cCHH=cCHH.3x.100.50p.all)
inNames <- list("fullCpG", "aCpG", "bCpG", "cCpG", "fullCHG", "aCHG", "bCHG", "hCHG", "fullCHH", "aCHH", "cCHH")

outRegion <- Map(function(diffRegions, diffNames) {
  
  print(diffNames)
  
  x <- as(diffRegions, "GRanges") # convert DM tile object to GRanges from methylKit
  y <- reduce(x, min.gapwidth=101) # join tiles, allowing one 100bp tile gap
  z <- y[width(y) > 300] # obtain only joined tiles of a minimum width
  
  subset1 <- subsetByOverlaps(z, up, invert=TRUE) # remove upstream genic flanks
  subset2 <- subsetByOverlaps(subset1, down,invert=TRUE) # remove downstream genic flanks
  subset3 <- subsetByOverlaps(subset2, geneObj,invert=TRUE) # remove genes
  subset4 <- subsetByOverlaps(subset3, te.file, invert=TRUE) # remove TEs from most recent annotation release
  subset5 <- subsetByOverlaps(subset4, legacy.file, invert=TRUE) # remove TEs from legacy file
  
  out <- resize(subset5, width=width(subset5)+1000, fix="center") # expand DMR 500bp in both directions
  
  outTable <- as.data.frame(out) # convert to dataframe
  outTable$seqnames <- as.numeric(as.character(outTable$seqnames)) # make chromosome numeric, not factor
  outTable$tempID <- rownames(outTable) # get unique number within comparison
  outTable$seqID <- paste0("> ",  outTable$tempID, " ", diffNames, " ", outTable$seqnames, ":", outTable$start, "-", outTable$end) # get fasta header line
  outTable$sequence <- as.character(subseq(refGenome[outTable$seqnames], start=outTable$start, end=outTable$end)) # get sequence data
  X <- data.frame(seqID=outTable$seqID, sequence=outTable$sequence) # set up dataframe
  D <- do.call(rbind, lapply(seq(nrow(X)), function(i) t(X[i, ]))) # fill dataframe
  write.csv(outTable, file=paste0("/home/nia/project/nia/results/2.6-novelRegions/csv-", diffNames, "-novel.csv"), row.names=F, quote=F) # export dataframe
  
  return(out)
  
}, inFiles, inNames)
