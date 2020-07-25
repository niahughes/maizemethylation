###################################################################
### TE ANNOTATION OF DMRs
### R script to annotate DM tiles with transposable elements
###################################################################

# set up environment
setwd("/scratch/nia/manuFinal")
load(file=".RData")
library(methylKit)
library(GenomicRanges)
library(genomation)
library(rtracklayer)

# Get annotation objects --------------------------------------------------
te.file <- import.gff3("/home/nia/project/nia/annotation/B73.structuralTEv2.fulllength.2018-09-19.gff3", format="gff3")
te.file <- keepStandardChromosomes(te.file, pruning.mode="coarse")
te.df <- as.data.frame(te.file)
te.df$code <- ""
te.df2 <- transform(te.df, code = substr(te.df$ID,1,3))
teObj <- makeGRangesFromDataFrame(te.df2, keep.extra.columns = TRUE)

legacy.file=keepStandardChromosomes(import.gff("/home/nia/project/nia/annotation/legacy.gff3", format="gff3"), pruning.mode="coarse")
mite = legacy.file[legacy.file$source=="detectMITE"]

# Function ----------------------------------------------------------------
annTEs <- function(diffMethTiles, compName) { # annotate with most recent annotation release
  x <- as(diffMethTiles, "GRanges") # convert to GRanges
  
  ann <- mergeByOverlaps(x, teObj) # merge with feature
  ann.df <- as.data.frame(ann) # convert to dataframe
  ann.df <- ann.df[ c(1:3,5,9:11,17,18,21) ] # remove extraneous columns
  names(ann.df) <- c("chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff", "source", "type", "ID", "Name", "code") # name columns
  
  write.table(file = paste0("/home/nia/project/nia/results/2.4-teAnnotate/TE-",compName,".csv"), 
              ann.df, sep = ",", row.names = F,quote=FALSE)
  
  return(ann.df)
}

annMITEs <- function(diffMethTiles, compName) { # annotate with MITEs from legacy annotation file
  x <- as(diffMethTiles, "GRanges") # convert to GRanges
  
  ann <- mergeByOverlaps(x, mite) # merge with feature
  ann.df <- as.data.frame(ann) # convert to dataframe
  ann.df <- transform(ann.df, code = substr(ann.df$ID,1,3))
  ann.df <- ann.df[ c(1:3,4,9:11,17,18,21,29) ] # remove extraneous columns
  names(ann.df) <- c("chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff", "source", "type", "ID","code") # name columns
  
  write.table(file = paste0("/home/nia/project/nia/results/2.4-teAnnotate/MITE-",compName,".csv"), 
              ann.df, sep = ",", row.names = F,quote=FALSE)
  
  return(ann.df)
}



# Run ---------------------------------------------------------------------
comps <- list("full.CpG", "A.CpG", "B.CpG", "C.CpG", 
              "full.CHG", "A.CHG", "B.CHG", "C.CHG", 
              "full.CHH", "A.CHH", "B.CHH", "C.CHH")
DM.list <- list(full.CpG=CpG.3x.100.oChi.50p.all, A.CpG=aCpG.3x.100.50p.all, B.CpG=bCpG.3x.100.50p.all, C.CpG=cCpG.3x.100.50p.all,
                full.CHG=CHG.3x.100.oChi.50p.all, A.CHG=aCHG.3x.100.50p.all, B.CHG=bCHG.3x.100.50p.all, C.CHG=hCHG.3x.100.50p.all,
                full.CHH=CHH.3x.100.oChi.50p.all, A.CHH=aCHH.3x.100.50p.all, B.CHH=bCHH.3x.100.50p.all, C.CHH=cCHH.3x.100.50p.all)
# annotate each set of differentially methylated tiles with TEs
teList <- Map(annTEs, diffMethTiles=DM.list, compName=comps)
miteList <- Map(annMITEs, diffMethTiles=DM.list, compName=comps)
