###################################################################
### GENIC ANNOTATION OF DMRs
### R script to annotate DM tiles with genes
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(methylKit)
library(GenomicRanges)
library(genomation)

# Get annotation object ---------------------------------------------------
gene.obj=readTranscriptFeatures("/home/nia/project/nia/annotation/geneAnnotationb73.bed", unique.prom=FALSE,  up.flank=2000, down.flank=0)
gene.obj$TSSes <- NULL
gene.types <- list("exons","introns","promoters")

# Function ----------------------------------------------------------------
annGenes <- function(diffMethTiles, compName) {
  x <- as(diffMethTiles, "GRanges") # convert to GRanges
  
  out <- lapply(gene.types, function(y) { # for each in exons, introns, promoters
    ann <- mergeByOverlaps(x, gene.obj[[y]]) # merge with feature
    ann.dedup <- ann[!duplicated(ann$x),] # deduplicate to remove doubled tiles
    ann.dedup$name <- gsub("transcript:|_.*$", "", ann.dedup$name) # remove "transcript" and ID from gene name
    ann.df <- as.data.frame(ann.dedup) # convert to dataframe
    ann.df <- ann.df[ -c(4,6:8,12:19) ] # remove extraneous/duplicated columns
    colnames(ann.df)<-c("chr", "start","end","strand","pvalue", "qvalue", "meth.diff", "name") # name columns
    ann.df$feature <- y # add column with feature type
    return(ann.df)
  })
  
  outList <- do.call(rbind, out) # combine all mapped features to one list
  rownames(outList) <- NULL # remove strange rownames
  
  write.table(file = paste0("/home/nia/project/nia/results/2.3-geneAnnotate/gene-", compName, ".csv"), 
              outList, sep = ",", row.names = F,quote=FALSE)
  
  return(outList)
}

# Run ---------------------------------------------------------------------
comps <- list("full.CpG", "A.CpG", "B.CpG", "C.CpG", 
              "full.CHG", "A.CHG", "B.CHG", "C.CHG", 
              "full.CHH", "A.CHH", "B.CHH", "C.CHH")
DM.list <- list(full.CpG=CpG.3x.100.oChi.50p.all, A.CpG=aCpG.3x.100.50p.all, B.CpG=bCpG.3x.100.50p.all, C.CpG=cCpG.3x.100.50p.all,
                full.CHG=CHG.3x.100.oChi.50p.all, A.CHG=aCHG.3x.100.50p.all, B.CHG=bCHG.3x.100.50p.all, C.CHG=hCHG.3x.100.50p.all,
                full.CHH=CHH.3x.100.oChi.50p.all, A.CHH=aCHH.3x.100.50p.all, B.CHH=bCHH.3x.100.50p.all, C.CHH=cCHH.3x.100.50p.all)
# annotate each set of differentially methylated tiles with genes
geneList <- Map(annGenes, diffMethTiles=DM.list, compName=comps)